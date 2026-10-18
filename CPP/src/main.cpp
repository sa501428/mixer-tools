#include <algorithm>

#include "slice/slice.hpp"
#include <iostream>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <thread>

namespace {
void usage() {
    std::cout << "Mixer Tools SLICe C++\n"
      "Usage: slice [slice] [options] <file.hic[,file2.hic]> <K0,KF,nK> <outfolder> <prefix>\n"
      "  K0 minimum K; KF exclusive upper bound; nK attempts per K\n"
      "  -r, --resolutions LIST   Resolution (default 100000; first value used)\n"
      "  -c, --chromosomes LIST   Selected autosomes\n"
      "  -z, --threads N          Workers (negative: hardware, zero: one)\n"
      "      --random-seeds LIST  Signed Java seeds; last value used\n"
      "      --subsample N        Java compatibility value (currently unused)\n"
      "      --encode-mode        Export all-K ENCODE BED\n"
      "  -v, --verbose            Verbose diagnostics\n"
      "  -h, --help               Help; -V, --version: source version\n"
      "C++ extensions exposing Java switches / dormant helpers:\n"
      "      --kmeans             Run k-means instead of default k-medians\n"
      "      --both               Run k-means followed by k-medians\n"
      "      --no-filter-outliers Disable global log-space filtering\n"
      "      --no-umap            Disable diagnostic projection\n"
      "      --has-translocation  Detect and mask suspected translocations\n"
      "  -k, --normalization NAME Override stage defaults (Java slice ignores -k)\n"
      "      --max-retries N      Bound Java's unbounded collapsed-cluster retries\n";
}
std::vector<std::string> split(const std::string& value) {
    std::vector<std::string> result;
    size_t start = 0;
    do {
        size_t end = value.find(',', start);
        auto part = value.substr(start, end == std::string::npos ? end : end - start);
        if (part.empty()) throw std::invalid_argument("Empty comma-separated value");
        result.push_back(part);
        if (end == std::string::npos) break;
        start = end + 1;
    } while (true);
    return result;
}
int64_t integer(const std::string& value) {
    size_t used = 0;
    auto n = std::stoll(value, &used);
    if (used != value.size()) throw std::invalid_argument("Invalid integer: " + value);
    return n;
}
int positive(const std::string& value) {
    auto n = integer(value);
    if (n < 1 || n > std::numeric_limits<int>::max()) throw std::invalid_argument("Expected positive 32-bit integer: " + value);
    return static_cast<int>(n);
}
}
int main(int argc, char** argv) {
    try {
        slice::Options options;
        options.threads = 1;
        std::vector<std::string> positional;
        bool endOptions = false;
        for (int i = 1; i < argc; ++i) {
            std::string arg(argv[i]);
            if (!endOptions && arg == "--") { endOptions = true; continue; }
            if (!endOptions && (arg == "-h" || arg == "--help")) { usage(); return 0; }
            if (!endOptions && (arg == "-V" || arg == "--version")) { std::cout << "Mixer Tools SLICe C++ 1.0\n"; return 0; }
            if (!endOptions && (arg == "-v" || arg == "--verbose")) { options.verbose = true; continue; }
            if (!endOptions && arg == "--encode-mode") { options.encode = true; continue; }
            if (!endOptions && arg == "--kmeans") { options.means = true; options.medians = false; continue; }
            if (!endOptions && arg == "--both") { options.means = options.medians = true; continue; }
            if (!endOptions && arg == "--no-filter-outliers") { options.filterOutliers = false; continue; }
            if (!endOptions && arg == "--no-umap") { options.project = false; continue; }
            if (!endOptions && arg == "--has-translocation") { options.detectTranslocations = true; continue; }
            if (!endOptions && arg.rfind('-', 0) == 0) {
                std::string name = arg.substr(0, arg.find('=')), value;
                const std::vector<std::string> known = {"-r", "--resolutions", "--resolution", "-c", "--chromosomes", "-z", "--threads", "--random-seeds", "--subsample", "-k", "--normalization", "--max-retries", "-m", "--matrix-window-width", "-w", "--window", "--seed", "--compare", "--type", "--corr"};
                if (std::find(known.begin(), known.end(), name) == known.end()) throw std::invalid_argument("Unknown option: " + name);
                if (arg.find('=') != std::string::npos) value = arg.substr(arg.find('=') + 1);
                else if (++i < argc) value = argv[i];
                else throw std::invalid_argument("Missing value for " + name);
                if (name == "--max-retries") options.maxRetries = positive(value);
                else if (name == "-r" || name == "--resolutions" || name == "--resolution") {
                    auto values = split(value); options.resolution = positive(values[0]);
                    if (values.size() > 1) std::cerr << "Only one resolution can be specified; using " << options.resolution << '\n';
                } else if (name == "-c" || name == "--chromosomes") options.chromosomes = split(value);
                else if (name == "-z" || name == "--threads") { auto n = integer(value); if (n > std::numeric_limits<int>::max()) throw std::invalid_argument("Thread count too large"); options.threads = n < 0 ? std::max(1U, std::thread::hardware_concurrency()) : std::max(1U, static_cast<unsigned>(n)); }
                else if (name == "-k" || name == "--normalization") {
                    if (value.empty()) throw std::invalid_argument("Normalization must not be empty");
                    options.normalization = value;
                }
                else if (name == "--random-seeds") for (const auto& seed : split(value)) options.seed = static_cast<uint64_t>(integer(seed));
                else if (name == "--subsample") { auto n = integer(value); if (n > std::numeric_limits<int>::max() || n < std::numeric_limits<int>::min()) throw std::invalid_argument("Subsample out of 32-bit range"); if (n > 0) options.subsample = static_cast<int>(n); }
                else { std::cerr << "Warning: " << name << " is unused by Java slice and has no effect\n"; }
            } else positional.push_back(arg);
        }
        if (!positional.empty() && positional[0] == "slice") positional.erase(positional.begin());
        if (positional.empty()) { usage(); return argc == 1 ? 0 : 1; }
        if (positional.size() != 4) throw std::invalid_argument("Expected file.hic, K0,KF,nK, output directory, and prefix (see --help)");
        options.files = split(positional[0]);
        auto k = split(positional[1]);
        if (k.size() != 3) throw std::invalid_argument("Cluster argument must be K0,KF,nK");
        options.firstK = positive(k[0]); options.endK = positive(k[1]); options.attempts = positive(k[2]);
        if (options.endK <= options.firstK) throw std::invalid_argument("KF must be greater than K0 (exclusive upper bound)");
        options.output = positional[2]; options.prefix = positional[3];
        if (options.prefix.find('/') != std::string::npos || options.prefix.find('\\') != std::string::npos) throw std::invalid_argument("Prefix must be a filename prefix");
        slice::run(options);
        return 0;
    } catch (const std::exception& e) { std::cerr << "slice: " << e.what() << '\n'; return 1; }
}
