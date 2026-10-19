#pragma once
#include "slice/slice.hpp"
#include "straw.h"
#include <fstream>
#include <functional>

// Internal native reader, chromosome compression and export operations.
namespace slice::detail {
struct Translocation {
    size_t a, b;
    int64_t startA, endA, startB, endB;
    double peakZ;
};
struct Input {
    std::string file, intra, inter, gw, coverage;
    StrawFileInfo info;
    std::vector<chromosome> chromosomes;
    std::vector<Translocation> translocations;
};
using Mask = std::vector<std::vector<bool>>;
struct Order { std::vector<int> assignments; size_t count; };
std::string canonical(std::string name);
std::string chooseNormLabel(const Input&,const std::vector<std::string>&);
size_t bins(const chromosome& chr, int resolution);
std::vector<Input> loadInputs(const Options& options);
void records(const Input& in, size_t a, size_t b, int resolution, const std::string& norm,
             const std::function<void(size_t, size_t, float)>& callback);
Mask badBins(const std::vector<Input>& inputs, int resolution, const Mask& excluded = {}, const std::filesystem::path& coverageReport = {});
Order orderChromosome(const Input& input, size_t a, const std::vector<bool>& bad,
                      const Options& options, Random& random);
CleanMatrix compressedTrans(const Input& input, const Mask& bad,
                            const std::vector<Order>& orders, int resolution, double minExposureFraction=.2, double memoryGB=4, bool legacy=false, unsigned threads=1, bool allowPartial=false);
std::vector<Translocation> findTranslocations(const Input& input, const Mask& bad, int resolution,
    int& detectionResolution, const std::filesystem::path& output, double memoryGB=4);
void saveBed(const std::filesystem::path&, const Input&, const CleanMatrix&, const std::vector<Cluster>&);
void saveNpy(const std::filesystem::path&, const std::vector<std::vector<double>>&);
void saveNpy(const std::filesystem::path&, const Matrix&);
void saveNpy(const std::filesystem::path&, const std::vector<std::vector<int32_t>>&);
std::ofstream outputFile(const std::filesystem::path& path, bool binary = false);
void bedLine(std::ostream& out, const chromosome& chr, int64_t start, int64_t end, int id);
}
