// Native diagnostic projection. Clustering never uses the UMAP coordinates.
#include "slice/slice.hpp"
#include "slice/parallel.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

namespace slice {
Matrix cosineFeatures(const Matrix& data, uint64_t seed, unsigned threads) {
    if (data.empty()) throw std::runtime_error("Empty projection input");
    Matrix centers;
    if (data.size() == 50) centers = data;
    else {
        Random quick(seed);
        auto groups = kmeans(data, std::min<int>(50, data.size()), 20, quick.nextLong(), threads, true);
        // QuickCentroids always returns MEAN profiles, even for k-medians groups.
        for (const auto& group : groups) if (group.members.size() > 5) {
            Row center(data[0].size(), 0); std::vector<size_t> counts(center.size(), 0);
            for (size_t i : group.members) for (size_t j = 0; j < center.size(); ++j)
                if (std::isfinite(data[i][j])) { center[j] += data[i][j]; ++counts[j]; }
            for (size_t j = 0; j < center.size(); ++j) center[j] = counts[j] ? center[j] / counts[j] : nan;
            centers.push_back(std::move(center));
        }
    }
    if (centers.empty()) throw std::runtime_error("No projection centroid has more than five rows; use --no-umap for small datasets");
    Matrix result(data.size(), Row(centers.size()));
    detail::parallelFor(data.size(), threads, [&](size_t i) {
        for (size_t c = 0; c < centers.size(); ++c) {
            double dot=0, xx=1e-100, yy=1e-100;
            for (size_t j = 0; j < centers[c].size(); ++j) {
                float product = data[i][j] * centers[c][j];
                if (!std::isnan(product)) { dot += product; xx += data[i][j]*data[i][j]; yy += centers[c][j]*centers[c][j]; }
            }
            float r = static_cast<float>(std::clamp(dot / std::sqrt(xx*yy), double(-.99f), double(.99f)));
            result[i][c] = static_cast<float>(std::log(static_cast<double>(1.f+r))-std::log(static_cast<double>(1.f-r)))/2.f;
        }
    });
    return result;
}

Matrix umap(const Matrix& data, unsigned threads) {
    // Exact Euclidean nearest neighbors, smooth-kNN local scaling, fuzzy union,
    // and negative-sampling cross-entropy optimization. Same n_neighbors=50,
    // min_dist=.5, spread=1, negative_sample_rate=5 and seed=0 as UmapProjection.
    // Native random initialization replaces tagbio's spectral/NN-descent code.
    const size_t n=data.size();
    if (n < 3) throw std::runtime_error("UMAP needs at least three points");
    if (data[0].size() <= 2) {
        Matrix answer=data; for (auto& row : answer) row.resize(2, 0); return answer;
    }
    const size_t k=std::min<size_t>(50,n-1);
    // UMAP counts the query point itself among its k nearest neighbors.
    const size_t partners=k-1;
    std::vector<std::vector<std::pair<size_t,double>>> local(n);
    detail::parallelFor(n, threads, [&](size_t i) {
        std::vector<std::pair<double,size_t>> nearest;
        nearest.reserve(n-1);
        for (size_t j=0;j<n;++j) if (i!=j) {
            double d=0;
            for (size_t c=0;c<data[i].size();++c) {
                if (!std::isfinite(data[i][c]) || !std::isfinite(data[j][c])) throw std::runtime_error("Nonfinite UMAP feature");
                double delta=data[i][c]-data[j][c]; d+=delta*delta;
            }
            nearest.emplace_back(std::sqrt(d),j);
        }
        std::partial_sort(nearest.begin(),nearest.begin()+partners,nearest.end());
        double rho=0, mean=0;
        for (size_t c=0;c<partners;++c) { mean+=nearest[c].first; if (!rho && nearest[c].first>0) rho=nearest[c].first; }
        double lo=0,hi=std::numeric_limits<double>::infinity(),sigma=1;
        for (int step=0;step<64;++step) {
            double sum=0; for (size_t c=0;c<partners;++c) sum+=std::exp(-std::max(0.,nearest[c].first-rho)/sigma);
            if (std::abs(sum-std::log2(double(k)))<1e-5) break;
            if (sum>std::log2(double(k))) { hi=sigma; sigma=(lo+hi)/2; }
            else { lo=sigma; sigma=std::isinf(hi)?sigma*2:(lo+hi)/2; }
        }
        sigma=std::max(sigma,1e-3*mean/k);
        for (size_t c=0;c<partners;++c) local[i].emplace_back(nearest[c].second,std::exp(-std::max(0.,nearest[c].first-rho)/std::max(sigma,1e-12)));
    });
    std::map<std::pair<size_t,size_t>, std::pair<double,double>> graph;
    for (size_t i=0;i<n;++i) for (auto [j,w]:local[i]) {
        auto& edge=graph[{std::min(i,j),std::max(i,j)}]; if (i<j) edge.first=w; else edge.second=w;
    }
    struct Edge { size_t i,j; double period,next,negativeNext; };
    std::vector<Edge> edges;
    const int epochs=n<=10000?500:200;
    for (auto [key,w]:graph) {
        double weight=w.first+w.second-w.first*w.second;
        if (weight>=1./epochs) for (int direction=0;direction<2;++direction)
            edges.push_back({direction?key.second:key.first,direction?key.first:key.second,1./weight,1./weight,1./weight/5});
    }
    Random random(0); Matrix embedding(n,Row(2));
    for (auto& row:embedding) for (float& x:row) x=float(random.nextInt(2000001))/100000.f-10.f;
    // Least-squares curve fit values for min_dist=.5, spread=1.
    constexpr double a=.5830300199,b=1.3341669931;
    for (int epoch=0;epoch<epochs;++epoch) {
        double alpha=1.-double(epoch)/epochs;
        for (auto& edge:edges) if (edge.next<=epoch) {
            auto& x=embedding[edge.i]; auto& y=embedding[edge.j];
            double d=0; for (int c=0;c<2;++c) d+=(x[c]-y[c])*(x[c]-y[c]);
            double coefficient=d>0?-2*a*b*std::pow(d,b-1)/(a*std::pow(d,b)+1):0;
            for (int c=0;c<2;++c) { float g=static_cast<float>(std::clamp(coefficient*(x[c]-y[c]),-4.,4.)*alpha); x[c]+=g; y[c]-=g; }
            edge.next+=edge.period;
            int negatives=static_cast<int>((epoch-edge.negativeNext)/(edge.period/5));
            for (int q=0;q<negatives;++q) {
                size_t j=static_cast<size_t>(random.nextInt(static_cast<int>(n))); if (j==edge.i) continue;
                double square=0; for (int c=0;c<2;++c) square+=(x[c]-embedding[j][c])*(x[c]-embedding[j][c]);
                double repulse=2*b/((.001+square)*(a*std::pow(square,b)+1));
                for (int c=0;c<2;++c) x[c]+=static_cast<float>(std::clamp(repulse*(x[c]-embedding[j][c]),-4.,4.)*alpha);
            }
            edge.negativeNext+=negatives*edge.period/5;
        }
    }
    return embedding;
}
}
