#pragma once
#include <cstdint>
#include <filesystem>
#include <limits>
#include <string>
#include <vector>

namespace slice {
using Row = std::vector<float>;
using Matrix = std::vector<Row>;
constexpr float nan = std::numeric_limits<float>::quiet_NaN();

// java.util.Random's 48-bit generator keeps the original seed hierarchy.
class Random {
    uint64_t state_;
    uint32_t next(unsigned bits);
public:
    explicit Random(uint64_t seed) : state_((seed ^ 0x5DEECE66DULL) & ((1ULL << 48) - 1)) {}
    uint64_t nextLong();
    int nextInt(int bound);
};
struct Options {
    std::vector<std::string> files, chromosomes;
    int resolution = 100000, firstK = 2, endK = 12, attempts = 3;
    unsigned threads = 1;
    uint64_t seed = 22871;
    bool verbose = false;
    std::filesystem::path output;
    std::string prefix;
    std::string normalization; // Empty selects the original per-stage defaults.
    bool encode = false, means = false, medians = true;
    bool filterOutliers = true, project = true, detectTranslocations = false;
    int subsample = 100, maxRetries = 100;

};
struct Interval {
    size_t chromosome;
    int64_t start, end;
};
struct Cluster {
    Row center;
    std::vector<size_t> members;
};
struct CleanMatrix {
    Matrix values;
    std::vector<int> weights;
    std::vector<Interval> intervals;
};
float correlation(const Row& x, const Row& y, bool arc = false);
double squaredDistance(const Row& x, const Row& y);
std::vector<Cluster> kmeans(const Matrix& data, int k, int iterations, uint64_t seed, unsigned threads, bool medians = false);
// SliceMatrixCleaner: global cleanup, sparse-row removal, z-score and weights.
void clean(CleanMatrix& data, bool filterOutliers = true);
void run(const Options& options);
double absoluteDistance(const Row& x, const Row& y);
void scaleWeights(CleanMatrix& data);
Matrix cosineFeatures(const Matrix& data, uint64_t seed, unsigned threads);
Matrix umap(const Matrix& data, unsigned threads);
}
