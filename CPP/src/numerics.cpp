// Port of the SLICE numerical routines. See LICENSE for the original MIT license.
#include "slice/slice.hpp"
#include "slice/parallel.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace slice {
uint32_t Random::next(unsigned bits) {
    state_ = (state_ * 0x5DEECE66DULL + 0xB) & ((1ULL << 48) - 1);
    return static_cast<uint32_t>(state_ >> (48 - bits));
}
uint64_t Random::nextLong() {
    const auto high = static_cast<int32_t>(next(32));
    const auto low = static_cast<int32_t>(next(32));
    return (static_cast<uint64_t>(high) << 32) + static_cast<uint64_t>(low);
}
int Random::nextInt(int bound) {
    if (bound <= 0) throw std::invalid_argument("Random bound must be positive");
    if ((bound & -bound) == bound) return static_cast<int>((static_cast<int64_t>(bound) * next(31)) >> 31);
    uint32_t bits, val;
    do { bits = next(31); val = bits % bound; } while (bits - val + (bound - 1U) >= 0x80000000U);
    return static_cast<int>(val);
}
float correlation(const Row& x, const Row& y, bool arc) {
    double sx = 0, sy = 0;
    size_t n = 0;
    for (size_t j = 0; j < x.size(); ++j) if (std::isfinite(x[j]) && std::isfinite(y[j])) {
        sx += x[j]; sy += y[j]; ++n;
    }
    if (!n) return nan;
    double dot = 0, xx = 0, yy = 0;
    for (size_t j = 0; j < x.size(); ++j) if (std::isfinite(x[j]) && std::isfinite(y[j])) {
        double a = x[j] - sx / n, b = y[j] - sy / n;
        dot += a * b; xx += a * a; yy += b * b;
    }
    if (xx == 0 || yy == 0) return nan;
    float r = static_cast<float>(dot / std::sqrt(xx * yy));
    if (arc) {
        r = std::clamp(r, -.99f, .99f);
        r = static_cast<float>(std::log(static_cast<double>(1.f + r)) - std::log(static_cast<double>(1.f - r))) / 2.f;
    }
    return r;
}
double squaredDistance(const Row& x, const Row& y) {
    double total = 0;
    size_t n = 0;
    for (size_t j = 0; j < x.size(); ++j) {
        float d = x[j] - y[j];
        if (std::isfinite(d)) { total += d * d; ++n; }
    }
    return n ? x.size() * total / n : std::numeric_limits<float>::max();
}

double absoluteDistance(const Row& x, const Row& y) {
    double sum = 0; size_t n = 0;
    for (size_t j = 0; j < x.size(); ++j) { float d = x[j] - y[j]; if (!std::isnan(d)) { sum += std::abs(d); ++n; } }
    return static_cast<float>((n ? sum / n : double(std::numeric_limits<float>::max())) * x.size());
}
void scaleWeights(CleanMatrix& data) {
    for (auto& row : data.values) for (size_t j = 0; j < row.size(); ++j) row[j] *= std::sqrt(data.weights[j]);
}
// Assignment is parallel; center reductions stay in row order so the result
// does not depend on scheduling or the requested thread count.
static void assign(const Matrix& data, const std::vector<Cluster>& clusters,
                   std::vector<int>& assignments, unsigned threads, bool medians) {
    detail::parallelFor(data.size(), threads, [&](size_t i) {
        double best = std::numeric_limits<double>::infinity();
        int chosen = -1;
        for (size_t c = 0; c < clusters.size(); ++c) if (!clusters[c].members.empty()) {
            double d = static_cast<float>(medians ? absoluteDistance(data[i], clusters[c].center) : std::sqrt(squaredDistance(data[i], clusters[c].center)));
            if (d < best) { best = d; chosen = static_cast<int>(c); }
        }
        assignments[i] = chosen;
    });
}
// Exact L1 lower bounds on each row's observed coordinates. If a centroid's
// support changes, its bounds are invalidated and the full distance is evaluated.
static void assignMedians(const Matrix& data,const std::vector<Cluster>& clusters,
                          std::vector<int>& assignments,unsigned threads,
                          std::vector<double>& lower,const std::vector<double>& movement,
                          const std::vector<size_t>& finite){
    size_t k=clusters.size();bool initial=lower.empty();if(initial)lower.assign(data.size()*k,0);
    detail::parallelFor(data.size(),threads,[&](size_t i){
        double factor=finite[i]?double(data[i].size())/finite[i]:1;
        if(!initial)for(size_t c=0;c<k;++c){double& bound=lower[i*k+c];
            bound=std::isfinite(bound)&&std::isfinite(movement[c])?std::max(0.,bound-movement[c]*factor-32*std::numeric_limits<float>::epsilon()*(1+bound)):0;}
        int chosen=assignments[i];double best=std::numeric_limits<double>::infinity();
        if(chosen>=0&&!clusters[size_t(chosen)].members.empty()){
            best=static_cast<float>(absoluteDistance(data[i],clusters[size_t(chosen)].center));lower[i*k+size_t(chosen)]=best;if(!std::isfinite(best))chosen=-1;}
        for(size_t c=0;c<k;++c)if(!clusters[c].members.empty()){
            if(!initial&&int(c)==chosen)continue;
            if(!initial&&lower[i*k+c]>best)continue;
            double distance=static_cast<float>(absoluteDistance(data[i],clusters[c].center));lower[i*k+c]=distance;
            if(distance<best||(std::isfinite(distance)&&distance==best&&(chosen<0||int(c)<chosen))){best=distance;chosen=int(c);}}
        assignments[i]=chosen;
    });
}
std::vector<Cluster> kmeans(const Matrix& data, int k, int iterations, uint64_t seed, unsigned threads, bool medians, bool bounded) {
    if (data.empty() || data[0].empty() || k < 1 || static_cast<size_t>(k) > data.size())
        throw std::invalid_argument("Not enough usable rows for requested cluster count");
    Random random(seed);
    std::vector<Cluster> clusters;
    Row closest(data.size(), std::numeric_limits<float>::max());
    size_t index = static_cast<size_t>(random.nextInt(static_cast<int>(data.size())));
    // K-medians repeatedly visits one coordinate across many members. A single
    // blocked transpose makes those scans contiguous and is shared by all centers.
    Matrix columns;
    if(medians){columns.assign(data[0].size(),Row(data.size()));
        for(size_t ib=0;ib<data.size();ib+=64)for(size_t jb=0;jb<data[0].size();jb+=64)
            for(size_t i=ib;i<std::min(data.size(),ib+64);++i)for(size_t j=jb;j<std::min(data[0].size(),jb+64);++j)columns[j][i]=data[i][j];}
    std::vector<int> assignments(data.size(), -1);
    for (int c = 0; c < k; ++c) {
        clusters.push_back({data[index], {index}});
        assignments[index] = c;
        detail::parallelFor(data.size(),threads,[&](size_t i){
            closest[i] = std::min(closest[i], static_cast<float>(medians ? absoluteDistance(data[i], data[index]) : std::sqrt(squaredDistance(data[i], data[index]))));});
        index = static_cast<size_t>(std::max_element(closest.begin(), closest.end()) - closest.begin());
    }
    std::vector<double> lower;std::vector<size_t> finite(data.size());
    if(medians&&bounded){for(size_t i=0;i<data.size();++i)for(float v:data[i])finite[i]+=std::isfinite(v);assignMedians(data,clusters,assignments,threads,lower,{},finite);}
    else assign(data, clusters, assignments, threads, medians);
    for (int iteration = 0; iteration < iterations; ++iteration) {
        for (auto& c : clusters) c.members.clear();
        for (size_t i = 0; i < data.size(); ++i) if (assignments[i] >= 0) clusters[assignments[i]].members.push_back(i);
        Matrix previousCenters;if(medians&&bounded)for(const auto& c:clusters)previousCenters.push_back(c.center);
        detail::parallelFor(clusters.size(),threads,[&](size_t ci){ auto& c=clusters[ci]; if (!c.members.empty()) {
            c.center.assign(data[0].size(), 0);
            if(!medians){
            std::vector<size_t> counts(c.center.size(), 0);
            for (auto i : c.members) for (size_t j = 0; j < c.center.size(); ++j) if (std::isfinite(data[i][j])) {
                c.center[j] += data[i][j]; ++counts[j];
            }
            for (size_t j = 0; j < c.center.size(); ++j) c.center[j] = counts[j] ? c.center[j] / counts[j] : nan;
            }
            Row values;if(medians)values.reserve(c.members.size());
            if (medians) for (size_t j = 0; j < c.center.size(); ++j) {
                values.clear();
                for (auto i : c.members) if (std::isfinite(columns[j][i])) values.push_back(columns[j][i]);
                if (values.empty()) { c.center[j] = nan; continue; }
                size_t mid = values.size() / 2;
                std::nth_element(values.begin(), values.begin()+mid, values.end());
                float value = values[mid];
                if (!(values.size() % 2)) value = (value + *std::max_element(values.begin(), values.begin()+mid)) / 2.f;
                c.center[j] = value;
            }
        }});
        auto before = assignments;
        if(medians&&bounded){std::vector<double> movement(clusters.size());
            for(size_t c=0;c<clusters.size();++c)for(size_t j=0;j<clusters[c].center.size();++j){double old=previousCenters[c][j],now=clusters[c].center[j];
                if(!std::isfinite(old)||!std::isfinite(now)){movement[c]=INFINITY;break;}
                movement[c]+=std::abs(old-now)+32*std::numeric_limits<float>::epsilon()*(std::abs(old)+std::abs(now)+1);}
            assignMedians(data,clusters,assignments,threads,lower,movement,finite);
        }else assign(data,clusters,assignments,threads,medians);
        if (before == assignments) break;
    }
    for (auto& c : clusters) c.members.clear();
    for (size_t i = 0; i < data.size(); ++i) if (assignments[i] >= 0) clusters[assignments[i]].members.push_back(i);
    clusters.erase(std::remove_if(clusters.begin(), clusters.end(), [](const Cluster& c) { return c.members.empty(); }), clusters.end());
    return clusters;
}

static std::pair<double, double> positiveStats(const Matrix& matrix) {
    double total = 0, squares = 0;
    size_t count = 0;
    for (const auto& row : matrix) for (float x : row) if (std::isfinite(x) && x > 0) { total += x; ++count; }
    if (!count) throw std::runtime_error("No positive inter-chromosomal values remain after filtering");
    double mean = total / count;
    for (const auto& row : matrix) for (float x : row) if (std::isfinite(x) && x > 0) squares += (x - mean) * (x - mean);
    return {mean, std::sqrt(squares / count)};
}
void clean(CleanMatrix& data, bool filterOutliers) {
    auto& m = data.values;
    if (m.empty() || m[0].empty() || data.weights.size() != m[0].size() || data.intervals.size() != m.size())
        throw std::invalid_argument("Empty or inconsistent SLICE matrix");
    if (filterOutliers) {
    for (auto& row : m) for (size_t j = 0; j < row.size(); ++j)
        row[j] = std::isfinite(row[j]) && row[j] >= 1e-20f && data.weights[j] > 0
               ? static_cast<float>(std::log(static_cast<double>(1.f + row[j] / data.weights[j]))) : nan;
    auto [mu, sd] = positiveStats(m);
    if (sd > 0) for (auto& row : m) for (float& x : row) if (std::isfinite(x) && x > 0 && (x - mu) / sd > 5) x = nan;
    auto stats = positiveStats(m); mu = stats.first; sd = stats.second;
    for (auto& row : m) for (float& x : row) if (std::isfinite(x)) {
        if (sd > 0 && std::abs((x - mu) / sd) > 2) x = nan;
        else x = static_cast<float>(std::expm1(static_cast<double>(x)));
    }
    }
    Matrix kept;
    std::vector<Interval> intervals;
    const auto limit = static_cast<size_t>(std::ceil(.7f * m[0].size()));
    for (size_t i = 0; i < m.size(); ++i) {
        size_t invalid = static_cast<size_t>(std::count_if(m[i].begin(), m[i].end(), [](float x) { return !std::isfinite(x); }));
        if (invalid <= limit) { kept.push_back(std::move(m[i])); intervals.push_back(data.intervals[i]); }
    }
    m = std::move(kept); data.intervals = std::move(intervals);
    if (m.empty()) throw std::runtime_error("All rows failed sparsity filtering; use a coarser resolution");
    for (size_t j = 0; j < m[0].size(); ++j) {
        double sum = 0, squares = 0; size_t n = 0;
        for (const auto& row : m) if (std::isfinite(row[j]) && row[j] > 0) { sum += row[j]; ++n; }
        float mean = static_cast<float>(sum / std::max<size_t>(1, n));
        for (const auto& row : m) if (std::isfinite(row[j]) && row[j] > 0) { float d = row[j] - mean; squares += d * d; }
        float stddev = static_cast<float>(std::sqrt(squares / std::max<size_t>(1, n)));
        for (auto& row : m) if (std::isfinite(row[j]))
            row[j] = stddev > 0 ? static_cast<float>((row[j] - mean) / stddev) : nan;
    }
    scaleWeights(data);
    for (const auto& row : m)
        if (std::none_of(row.begin(),row.end(),[](float value) { return std::isfinite(value); }))
            throw std::runtime_error("A genomic row has no variable features after z-scoring; use a coarser resolution");
}
}
