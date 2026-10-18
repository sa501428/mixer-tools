// Native SLICE pipeline; all .hic access goes through straw C++.
#include "slice/slice.hpp"
#include "slice/parallel.hpp"
#include "pipeline_internal.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <cctype>

namespace slice {
namespace detail {
std::string canonical(std::string name) {
    for (char& c : name) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if (name.rfind("CHR", 0) == 0) name.erase(0, 3);
    return name;
}
size_t bins(const chromosome& chr, int resolution) {
    if (chr.length <= 0) throw std::runtime_error("Invalid chromosome length: " + chr.name);
    const uint64_t count = static_cast<uint64_t>(chr.length / resolution) + 1;
    if (count > static_cast<uint64_t>(std::numeric_limits<int>::max())) throw std::runtime_error("Too many chromosome bins");
    return static_cast<size_t>(count);
}
void requireResolution(const Input& input, int resolution) {
    if (std::find(input.info.bpResolutions.begin(), input.info.bpResolutions.end(), resolution) == input.info.bpResolutions.end())
        throw std::runtime_error(input.file + ": missing BP resolution " + std::to_string(resolution));
}
bool normalizationVector(const Input& input, const chromosome& chr, int resolution,
                         const std::string& norm, std::vector<double>& values) {
    // The local straw v9 reader throws for this missing capability, whereas
    // its v8/v10 paths return false. Preserve actual I/O/corruption failures.
    try {
        return getNormalizationVectorForFile(input.file, chr.name, resolution, norm, values);
    } catch (const std::runtime_error& error) {
        if (std::string(error.what()) == "v9 normalization vectors not found for requested chromosomes, unit and resolution") return false;
        throw;
    }
}
std::string chooseNormLabel(const Input& input, const std::vector<std::string>& candidates) {
    for (const auto& norm : candidates)
        if (norm == "NONE" || std::find(input.info.normalizations.begin(),input.info.normalizations.end(),norm)!=input.info.normalizations.end()) return norm;
    std::string names;
    for (const auto& candidate : candidates) { if (!names.empty()) names += "/"; names += candidate; }
    throw std::runtime_error(input.file + ": no " + names + " normalization in dataset");
}
std::string chooseNorm(const Input& input, int resolution, const std::vector<std::string>& candidates) {
    // NormalizationPicker chooses the first advertised label, not the first
    // label with usable vectors at this resolution. Missing vectors are errors.
    auto norm=chooseNormLabel(input,candidates);
    for (const auto& chr : input.chromosomes) {
        std::vector<double> vector;
        if (!normalizationVector(input,chr,resolution,norm,vector) || vector.empty())
            throw std::runtime_error(input.file+": missing "+norm+" normalization vector for "+chr.name+" at "+std::to_string(resolution));
    }
    return norm;
}
// Prepared queries validate the matrix and normalization and return normalized
// counts in base-pair coordinates. Convert them once at the reader boundary.
void records(const Input& in, size_t a, size_t b, int resolution, const std::string& norm,
             const std::function<void(size_t, size_t, float)>& callback) {
    try {
        StrawPreparedQuery query(in.file, "observed", norm, in.chromosomes[a].name, in.chromosomes[b].name, "BP", resolution);
        size_t count = 0;
        query.streamWindow(0, in.chromosomes[a].length, 0, in.chromosomes[b].length, [&](const contactRecord& r) {
            if (r.binX < 0 || r.binY < 0) throw std::runtime_error("Negative contact coordinate");
            size_t x = static_cast<size_t>(r.binX / resolution), y = static_cast<size_t>(r.binY / resolution);
            if (x >= bins(in.chromosomes[a], resolution) || y >= bins(in.chromosomes[b], resolution))
                throw std::runtime_error("Contact coordinate outside chromosome");
            ++count;
            if (a != b) for (const auto& region : in.translocations) {
                if (region.a == a && region.b == b && r.binX >= region.startA && r.binX < region.endA &&
                    r.binY >= region.startB && r.binY < region.endB) return;
            }
            callback(x, y, r.counts);
        });
        if (!count) throw std::runtime_error("Matrix contains no contact records");
    } catch (const std::exception& e) {
        throw std::runtime_error(in.file + ": " + in.chromosomes[a].name + " x " + in.chromosomes[b].name +
                                 " at " + std::to_string(resolution) + " (" + norm + "): " + e.what());
    }
}
std::pair<float, float> coverageStats(const Matrix& data) {
    double sum = 0, square = 0; size_t n = 0;
    for (const auto& row : data) for (float x : row) if (x > 1e-10f && x < .9f * std::numeric_limits<float>::max()) { sum += x; ++n; }
    if (!n) throw std::runtime_error("No inter-chromosomal coverage; use a coarser resolution");
    float mean = static_cast<float>(sum / n);
    for (const auto& row : data) for (float x : row) if (x > 1e-10f && x < .9f * std::numeric_limits<float>::max()) {
        float diff = x - mean; square += diff * diff;
    }
    return {mean, static_cast<float>(std::sqrt(square / n))};
}
Mask badBins(const std::vector<Input>& inputs, int resolution, const Mask& excluded) {
    Matrix sums, nonzeros;
    std::vector<std::vector<size_t>> lengths;
    for (const auto& chr : inputs[0].chromosomes) {
        sums.emplace_back(bins(chr, resolution), 0); nonzeros.emplace_back(bins(chr, resolution), 0);
        lengths.emplace_back(bins(chr, resolution), 0);
    }
    for (const auto& input : inputs) for (size_t a = 0; a < sums.size(); ++a) for (size_t b = a + 1; b < sums.size(); ++b) {
        records(input, a, b, resolution, input.coverage.empty() ? input.intra : input.coverage, [&](size_t x, size_t y, float count) {
            if (!excluded.empty() && (excluded[a][x] || excluded[b][y])) return;
            float value = static_cast<float>(std::log(static_cast<double>(count + 1.f)));
            if (!std::isfinite(value) || value < 1e-10f) return;
            sums[a][x] += value; sums[b][y] += value; ++nonzeros[a][x]; ++nonzeros[b][y];
        });
        size_t validA = excluded.empty() ? sums[a].size() : std::count(excluded[a].begin(), excluded[a].end(), false);
        size_t validB = excluded.empty() ? sums[b].size() : std::count(excluded[b].begin(), excluded[b].end(), false);
        for (auto& length : lengths[a]) length += validB;
        for (auto& length : lengths[b]) length += validA;
        for (int direction = 0; direction < 2; ++direction) {
            size_t target = direction ? b : a, partner = direction ? a : b;
            std::map<size_t, std::vector<bool>> masks;
            for (const auto& region : input.translocations) if (region.a == a && region.b == b) {
                int64_t x0 = direction ? region.startB : region.startA, x1 = direction ? region.endB : region.endA;
                int64_t y0 = direction ? region.startA : region.startB, y1 = direction ? region.endA : region.endB;
                for (size_t x = static_cast<size_t>(x0 / resolution); x < sums[target].size() && int64_t(x) * resolution < x1; ++x) {
                    auto& mask = masks[x];
                    if (mask.empty()) mask.assign(sums[partner].size(), false);
                    for (size_t y = static_cast<size_t>(y0 / resolution); y < mask.size() && int64_t(y) * resolution < y1; ++y) mask[y] = true;
                }
            }
            for (const auto& row : masks) for (size_t y = 0; y < row.second.size(); ++y)
                if (row.second[y] && (excluded.empty() || !excluded[partner][y])) --lengths[target][row.first];
        }
    }
    for (size_t a = 0; a < sums.size(); ++a) for (size_t x = 0; x < sums[a].size(); ++x) {
        if (lengths[a][x]) { sums[a][x] /= lengths[a][x]; nonzeros[a][x] /= lengths[a][x]; }
        else { sums[a][x] = 0; nonzeros[a][x] = 0; }
    }
    auto [mean, sd] = coverageStats(sums);
    auto [nzmean, nzsd] = coverageStats(nonzeros);
    Mask bad;
    for (size_t a = 0; a < sums.size(); ++a) {
        bad.emplace_back(sums[a].size(), false);
        for (size_t x = 0; x < sums[a].size(); ++x)
            bad[a][x] = !(sums[a][x] > 0) || !(nonzeros[a][x] > 0) ||
                        (sd > 0 && (sums[a][x] - mean) / sd > 5) ||
                        (nzsd > 0 && (nonzeros[a][x] - nzmean) / nzsd < -3);
    }
    return bad;
}
Order orderChromosome(const Input& input, size_t a, const std::vector<bool>& bad,
                      const Options& options, Random& random) {
    int lowres = std::max(options.resolution, 100000), factor = lowres / options.resolution;
    size_t n = bins(input.chromosomes[a], lowres);
    Matrix oe(n, Row(n, 0));
    std::vector<double> expected;
    if (!getExpectedVectorForFile(input.file, input.chromosomes[a].name, lowres, input.intra, expected) || expected.empty())
        throw std::runtime_error("Missing " + input.intra + " expected vector for " + input.chromosomes[a].name + " at " + std::to_string(lowres));
    // ExpectedValueFunction.getCorrectedVersion: 5 Mb rolling median.
    // Reproduce QuickMedian's growing prefix and sliding-tail indexing.
    const size_t window = 5000000 / lowres;
    if (window > 0 && window < expected.size()) {
        if (2 * window >= expected.size()) throw std::runtime_error("Expected vector too short for Java's correction window");
        std::vector<double> smooth(expected.size()), values{expected[0]};
        auto median = [](std::vector<double> v) {
            size_t mid = v.size()/2; std::nth_element(v.begin(), v.begin()+mid, v.end());
            return v.size()%2 ? v[mid] : (v[mid]+*std::max_element(v.begin(),v.begin()+mid))/2;
        };
        size_t next = 0;
        for (size_t z=0; z<window; ++z) {
            smooth[z]=median(values); values.push_back(expected[next+1]); values.push_back(expected[next+2]); next+=2;
        }
        for (size_t z=window; z<expected.size(); ++z) {
            smooth[z]=median(values); values.erase(values.begin());
            if (window+1+z < expected.size()) values.push_back(expected[window+1+z]);
        }
        expected=std::move(smooth);
    }
    records(input, a, a, lowres, input.intra, [&](size_t x, size_t y, float observed) {
        size_t d = x > y ? x - y : y - x;
        double value = (observed + 1.) / (expected[std::min(d, expected.size() - 1)] + 1.);
        float v = std::isfinite(value) ? static_cast<float>(value) : 0.f;
        oe[x][y] = oe[y][x] = v;
    });
    std::vector<int> assignments(n, -5);
    for (size_t x = 0; x < bad.size(); ++x) if (bad[x]) assignments[x / factor] = -1;
    for (size_t x = 0; x < n; ++x) if (assignments[x] == -1) {
        std::fill(oe[x].begin(), oe[x].end(), nan);
        for (auto& row : oe) row[x] = nan;
    }
    // IntraMatrixCleaner.oeClean: exclude first three diagonals, zeros and O/E > 30.
    for (size_t x = 0; x < n; ++x) for (size_t y = 0; y < n; ++y)
        if ((x > y ? x-y : y-x) < 3 || std::abs(oe[x][y]) < 1e-10f || std::abs(oe[x][y]) > 30) oe[x][y] = nan;
    if (options.verbose) { std::vector<std::vector<double>> dump; for (auto row:oe) dump.emplace_back(row.begin(),row.end()); saveNpy(options.output/("cis_oe_"+input.chromosomes[a].name+".npy"),dump); }
    Matrix similarity(n, Row(n, nan));
    detail::parallelFor(n, options.threads, [&](size_t x) {
        if (assignments[x] < -2) {
            for (size_t y = x + 1; y < n; ++y) if (assignments[y] < -2) similarity[x][y] = similarity[y][x] = correlation(oe[x], oe[y], true);
        }
    });
    for (size_t x = 0; x < n; ++x) for (size_t y = 0; y < n; ++y)
        if ((x > y ? x-y : y-x) < static_cast<size_t>(5000000 / lowres)) similarity[x][y] = nan;
    // Release the dense O/E matrix before allocating k-means working data.
    Matrix{}.swap(oe);
    if (options.verbose) { std::vector<std::vector<double>> dump; for (auto row:similarity) dump.emplace_back(row.begin(),row.end()); saveNpy(options.output/("cis_corr_"+input.chromosomes[a].name+".npy"),dump); }
    Matrix usable;
    for (size_t x = 0; x < n; ++x) if (assignments[x] < -2) {
        auto invalid = std::count_if(similarity[x].begin(), similarity[x].end(), [](float v) { return !std::isfinite(v); });
        if (static_cast<float>(invalid) / n < .7) usable.push_back(similarity[x]);
    }
    int potential = std::max<int>(static_cast<int>(input.chromosomes[a].length / 50000000) + 7, 7);
    if (usable.size() < static_cast<size_t>(potential))
        throw std::runtime_error(input.chromosomes[a].name + ": too few usable intra-chromosomal profiles for compression");
    auto quickSeed = random.nextLong();
    if (options.verbose) std::cout << "Compression seed " << static_cast<int64_t>(quickSeed) << " K " << potential << '\n';
    Random quick(quickSeed);
    auto clusters = kmeans(usable, potential, 20, quick.nextLong(), options.threads);
    Matrix centroids;
    for (const auto& c : clusters) if (c.members.size() > 2) {
        Row center(n, 0); std::vector<size_t> counts(n, 0);
        for (size_t i : c.members) for (size_t j = 0; j < n; ++j) if (std::isfinite(usable[i][j])) {
            center[j] += usable[i][j]; ++counts[j];
        }
        for (size_t j = 0; j < n; ++j) center[j] = counts[j] ? center[j] / counts[j] : nan;
        centroids.push_back(std::move(center));
    }
    if (options.verbose) { std::vector<std::vector<double>> dump; for (auto row:centroids) dump.emplace_back(row.begin(),row.end()); saveNpy(options.output/("cis_centroids_"+input.chromosomes[a].name+".npy"),dump); }
    if (centroids.empty()) throw std::runtime_error(input.chromosomes[a].name + ": no compression centroid has more than two bins");
    std::vector<int> clustered(n, -1);
    detail::parallelFor(n, options.threads, [&](size_t x) {
        if (assignments[x] >= -2) return;
        float best = .2f;
        for (size_t c = 0; c < centroids.size(); ++c) {
            float r = correlation(centroids[c], similarity[x]);
            if (r > best) { best = r; clustered[x] = static_cast<int>(c); }
        }
    });
    for (size_t x = 0; x < n; ++x) if (assignments[x] < -2) assignments[x] = clustered[x];
    Order result{{}, centroids.size()};
    for (size_t x = 0; x < bad.size(); ++x) result.assignments.push_back(assignments[x / factor]);
    return result;
}
std::ofstream outputFile(const std::filesystem::path& path, bool binary) {
    std::ofstream out(path, binary ? std::ios::binary : std::ios::out);
    if (!out) throw std::runtime_error("Cannot write " + path.string());
    out.exceptions(std::ios::badbit | std::ios::failbit);
    return out;
}
const char* color(int id) {
    static const char* colors[] = {"0,0,0", "34,139,34", "152,251,152", "220,20,60", "255,255,0", "112,128,144", "75,0,130", "255,255,0", "255,0,0", "0,234,255", "170,0,255", "255,127,0", "191,255,0", "0,149,255", "255,0,170", "255,212,0", "106,255,0", "0,64,255", "237,185,185", "185,215,237", "231,233,185", "220,185,237", "185,237,224", "143,35,35", "35,98,143", "143,106,35", "107,35,143", "79,143,35", "115,115,115", "204,204,204"};
    return colors[std::abs(id) % 30];
}
void bedLine(std::ostream& out, const chromosome& chr, int64_t start, int64_t end, int id) {
    std::string name = chr.name;
    if (name.rfind("chr", 0) != 0) name = "chr" + name;
    out << name << '\t' << start << '\t' << end << '\t' << id << '\t' << id << "\t.\t" << start << '\t' << end << '\t' << color(id) << '\n';
}
void saveBed(const std::filesystem::path& path, const Input& input, const CleanMatrix& data, const std::vector<Cluster>& clusters) {
    std::vector<int> labels(data.values.size(), 0);
    for (size_t c = 0; c < clusters.size(); ++c) for (size_t i : clusters[c].members) labels[i] = static_cast<int>(c + 1);
    auto out = outputFile(path);
    for (size_t i = 0; i < labels.size();) {
        size_t end = i + 1;
        while (end < labels.size() && labels[end] == labels[i] && data.intervals[end].chromosome == data.intervals[i].chromosome && data.intervals[end].start == data.intervals[end - 1].end) ++end;
        bedLine(out, input.chromosomes[data.intervals[i].chromosome], data.intervals[i].start, data.intervals[end - 1].end, labels[i]);
        i = end;
    }
}
template<class T> void writeNpy(const std::filesystem::path& path, const std::vector<std::vector<T>>& values, const std::string& dtype) {
    if (values.empty() || values[0].empty()) throw std::runtime_error("Cannot write empty NPY matrix");
    auto out = outputFile(path, true);
    std::string header = "{'descr': '" + dtype + "', 'fortran_order': False, 'shape': (" + std::to_string(values.size()) + ", " + std::to_string(values[0].size()) + "), }";
    header.append((64 - ((10 + header.size() + 1) % 64)) % 64, ' '); header += '\n';
    if (header.size()>65535) throw std::runtime_error("NPY header too large");
    out.write("\x93NUMPY\x01\x00", 8);
    char size[2] = {static_cast<char>(header.size() & 255), static_cast<char>((header.size() >> 8) & 255)};
    out.write(size, 2); out.write(header.data(), static_cast<std::streamsize>(header.size()));
    const uint16_t endian = 1;
    for (const auto& row : values) {
        if (row.size()!=values[0].size()) throw std::runtime_error("Ragged NPY matrix");
        for (T x : row) {
            const char* bytes = reinterpret_cast<const char*>(&x);
            if (*reinterpret_cast<const char*>(&endian) == 1) out.write(bytes, sizeof x);
            else for (int i = sizeof(T)-1; i >= 0; --i) out.put(bytes[i]);
        }
    }
}
void saveNpy(const std::filesystem::path& path, const std::vector<std::vector<double>>& values) { writeNpy(path,values,"<f8"); }
void saveNpy(const std::filesystem::path& path, const Matrix& values) { writeNpy(path,values,"<f4"); }
void saveNpy(const std::filesystem::path& path, const std::vector<std::vector<int32_t>>& values) { writeNpy(path,values,"<i4"); }
std::vector<Input> loadInputs(const Options& options) {
    if (options.resolution <= 0 || options.firstK < 1 || options.endK <= options.firstK || options.attempts < 1 || !options.threads || options.files.empty())
        throw std::invalid_argument("Invalid SLICE options");
    const int lowres = std::max(options.resolution, 100000);
    if (lowres % options.resolution) throw std::invalid_argument("Resolution must divide 100000 when finer than 100000 BP");
    std::vector<Input> inputs;
    std::set<std::string> requested;
    for (const auto& name : options.chromosomes) requested.insert(canonical(name));
    for (const auto& file : options.files) {
        Input in; in.file = file; in.info = getFileInfo(file);
        if (inputs.empty()) {
            for (const auto& chr : in.info.chromosomes) {
                auto key = canonical(chr.name);
                if (key == "ALL" || key == "X" || key == "Y" || key == "M" || key == "MT") continue;
                if (requested.empty() || requested.count(key)) in.chromosomes.push_back(chr);
            }
            for (const auto& key : requested)
                if (std::none_of(in.chromosomes.begin(), in.chromosomes.end(), [&](const chromosome& c) { return canonical(c.name) == key; }))
                    throw std::runtime_error("Requested autosome not found: " + key);
            if (in.chromosomes.size() < 2) throw std::runtime_error("SLICE requires at least two autosomes with inter-chromosomal contacts");
        } else {
            if (in.info.genome != inputs[0].info.genome) throw std::runtime_error("Input genome identifiers do not match");
            for (const auto& reference : inputs[0].chromosomes) {
                auto found = std::find_if(in.info.chromosomes.begin(), in.info.chromosomes.end(), [&](const chromosome& c) { return canonical(c.name) == canonical(reference.name); });
                if (found == in.info.chromosomes.end() || found->length != reference.length) throw std::runtime_error("Input chromosomes or lengths do not match: " + reference.name);
                in.chromosomes.push_back(*found);
            }
        }
        requireResolution(in, options.resolution);
        requireResolution(in, lowres);
        in.intra = chooseNorm(in, options.resolution, options.normalization.empty()
            ? std::vector<std::string>{"SCALE", "KR"} : std::vector<std::string>{options.normalization});
        in.gw = chooseNormLabel(in, options.normalization.empty()
            ? std::vector<std::string>{"GW_SCALE", "GW_KR"} : std::vector<std::string>{options.normalization});
        in.coverage = in.intra;
        {
            // The same intra normalization must work for both coverage and O/E.
            for (const auto& chr : in.chromosomes) {
                std::vector<double> vec;
                if (!normalizationVector(in, chr, lowres, in.intra, vec) || vec.empty())
                    throw std::runtime_error("Missing " + in.intra + " normalization at compression resolution " + std::to_string(lowres));
            }
            in.inter = chooseNorm(in, options.resolution, options.normalization.empty()
                ? std::vector<std::string>{"INTER_SCALE", "INTER_KR"} : std::vector<std::string>{options.normalization});
        }
        if (options.verbose) std::cout << file << ": intra=" << in.intra << " inter=" << in.inter << " gw=" << in.gw << '\n';
        inputs.push_back(std::move(in));
    }
    return inputs;
}
CleanMatrix compressedTrans(const Input& input, const Mask& bad, const std::vector<Order>& orders, int resolution) {
    CleanMatrix data;
    std::vector<std::vector<int64_t>> rows;
    std::vector<size_t> offsets;
    size_t columns = 0;
    for (size_t a = 0; a < orders.size(); ++a) {
        offsets.push_back(columns); columns += orders[a].count;
        rows.emplace_back(bad[a].size(), -1);
        for (size_t x = 0; x < bad[a].size(); ++x) if (!bad[a][x]) {
            rows[a][x] = static_cast<int64_t>(data.intervals.size());
            int64_t start = static_cast<int64_t>(x) * resolution;
            data.intervals.push_back({a, start, std::min(start + resolution, input.chromosomes[a].length)});
        }
    }
    if (!columns || data.intervals.empty()) throw std::runtime_error("No matrix created; map is too sparse");
    data.weights.assign(columns, 0);
    data.values.assign(data.intervals.size(), Row(columns, 0));
    for (size_t a = 0; a < orders.size(); ++a) {
        for (size_t x = 0; x < bad[a].size(); ++x) if (!bad[a][x] && orders[a].assignments[x] >= 0)
            ++data.weights[offsets[a] + orders[a].assignments[x]];
        for (auto row : rows[a]) if (row >= 0)
            std::fill(data.values[row].begin() + offsets[a], data.values[row].begin() + offsets[a] + orders[a].count, nan);
    }
    for (size_t a = 0; a < orders.size(); ++a) for (size_t b = a + 1; b < orders.size(); ++b)
        records(input, a, b, resolution, input.inter, [&](size_t x, size_t y, float v) {
            if (!std::isfinite(v)) return;
            int ca = orders[a].assignments[x], cb = orders[b].assignments[y];
            if (rows[a][x] >= 0 && !bad[b][y] && cb >= 0) data.values[rows[a][x]][offsets[b] + cb] += v;
            if (rows[b][y] >= 0 && !bad[a][x] && ca >= 0) data.values[rows[b][y]][offsets[a] + ca] += v;
        });
    if (!input.translocations.empty()) {
        // Union excluded partner bins before counting, so overlapping regions
        // cannot subtract the same bin twice. Only affected target rows are built.
        using PartnerMask = std::map<size_t, std::vector<bool>>;
        std::map<std::pair<size_t, size_t>, PartnerMask> masks;
        for (const auto& region : input.translocations) for (int direction = 0; direction < 2; ++direction) {
            size_t a = direction ? region.b : region.a, b = direction ? region.a : region.b;
            int64_t x0 = direction ? region.startB : region.startA, x1 = direction ? region.endB : region.endA;
            int64_t y0 = direction ? region.startA : region.startB, y1 = direction ? region.endA : region.endB;
            for (size_t x = static_cast<size_t>(x0 / resolution); x < bad[a].size() && int64_t(x) * resolution < x1; ++x) {
                if (rows[a][x] < 0) continue;
                auto& mask = masks[{a,b}][x];
                if (mask.empty()) mask.assign(bad[b].size(), false);
                for (size_t y = static_cast<size_t>(y0 / resolution); y < mask.size() && int64_t(y) * resolution < y1; ++y) mask[y] = true;
            }
        }
        for (const auto& pair : masks) {
            const size_t a = pair.first.first, b = pair.first.second;
            for (const auto& row : pair.second) {
                std::vector<int> removed(orders[b].count, 0);
                for (size_t y = 0; y < row.second.size(); ++y)
                    if (row.second[y] && !bad[b][y] && orders[b].assignments[y] >= 0) ++removed[orders[b].assignments[y]];
                for (size_t g = 0; g < removed.size(); ++g) if (removed[g]) {
                    size_t j = offsets[b] + g; int eligible = data.weights[j] - removed[g];
                    auto& value = data.values[rows[a][row.first]][j];
                    value = eligible > 0 ? static_cast<float>(double(value) * data.weights[j] / eligible) : nan;
                }
            }
        }
    }
    return data;
}
}
}
