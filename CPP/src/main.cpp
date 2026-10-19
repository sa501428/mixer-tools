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
      "Usage: slice [options] FILE.hic OUTPUT [PREFIX] (automatic K=1..10)\n"
      "   or: slice [slice] [options] <file.hic[,file2.hic]> <K0,KF,nK> <outfolder> <prefix>\n"
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
      "      --umap               Enable exact diagnostic projection (quadratic; default off)\n"
      "      --no-umap            Disable diagnostic projection\n"
      "      --match-overlap F      Required mutual overlap in (0.5,1] (0.8)\n"
      "      --cis-view-weight F    Cis vs inter consensus weight (1)\n"
      "      --normalization-views inter|dual (inter)\n"
      "      --cluster-policy fixed|adaptive (fixed)\n"
      "      --target-clusters K  Target for adaptive search\n"
      "      --k-range MIN,MAX    Inclusive adaptive search bounds\n"
      "      --k-search screened|exhaustive Restart budget (default screened)\n"
    "      --k-selection bic|aic  Automatic selection criterion (default bic)\n"
      "      --max-iterations N    Per-fit bound (200)\n"
      "      --model kmeans|gmm  Adaptive continuous backend\n"
      "      --has-translocation  Detect and mask suspected translocations\n"
      "  -k, --normalization NAME Override stage defaults (Java slice ignores -k)\n"
      "      --cis-normalization NAME Cis/coverage norm (SCALE, fallback KR)\n"
      "      --trans-normalization NAME NONE (computed balance), AUTO, or stored norm\n"
      "      --post-balance       Explicitly rebalance stored-normalized trans inputs\n"
      "      --legacy-scale       Preserve cis SCALE throughout and legacy cleaning\n"
      "      --java-compatible    Stored stage norms and original Java cleaning\n"
      "      --row-standardization none|weighted-zscore (none)\n"
      "      --min-exposure-fraction F Minimum eligible target fraction (0.2)\n"
      "      --max-memory-gb F    Matrix workspace budget (4 GiB)\n"
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
        std::cout<<std::unitbuf;
        slice::Options options;
        options.threads = std::min(8u,std::max(1u,std::thread::hardware_concurrency()));
        std::vector<std::string> positional;
        bool endOptions = false;
        for (int i = 1; i < argc; ++i) {
            std::string arg(argv[i]);
            if (!endOptions && arg == "--") { endOptions = true; continue; }
            if (!endOptions && (arg == "-h" || arg == "--help")) { usage(); return 0; }
            if (!endOptions && (arg == "-V" || arg == "--version")) { std::cout << "Mixer Tools SLICe C++ 1.0\n"; return 0; }
            if (!endOptions && (arg == "-v" || arg == "--verbose")) { options.verbose = true; continue; }
            if (!endOptions && arg == "--post-balance") {options.forceBalance=true;continue;}
            if (!endOptions && arg == "--legacy-scale") { options.normalization="SCALE"; options.transNormalization="SCALE"; options.balanceTrans=false; options.legacyCleaning=true; continue; }
            if (!endOptions && arg == "--java-compatible") { options.javaCompatible=true; options.project=true; options.balanceTrans=false; options.legacyCleaning=true; continue; }
            if (!endOptions && arg == "--encode-mode") { options.encode = true; continue; }
            if (!endOptions && arg == "--kmeans") { options.means = true; options.medians = false; continue; }
            if (!endOptions && arg == "--both") { options.means = options.medians = true; continue; }
            if (!endOptions && arg == "--no-filter-outliers") { options.filterOutliers = false; continue; }
            if (!endOptions && arg == "--serial-attempts") {options.parallelAttempts=false;continue;}
            if (!endOptions && arg == "--no-distance-bounds") {options.distanceBounds=false;continue;}
            if (!endOptions && arg == "--umap") { options.project=true; continue; }
            if (!endOptions && arg == "--no-umap") { options.project = false; continue; }
            if (!endOptions && arg == "--has-translocation") { options.detectTranslocations = true; continue; }
            if (!endOptions && arg.rfind('-', 0) == 0) {
                std::string name = arg.substr(0, arg.find('=')), value;
                const std::vector<std::string> known = {"-r", "--resolutions", "--resolution", "-c", "--chromosomes", "-z", "--threads", "--random-seeds", "--subsample", "-k", "--normalization", "--max-retries", "--max-iterations", "--normalization-views", "--cluster-policy", "--model", "--k-range", "--target-clusters", "--k-selection", "--k-search", "--cis-view-weight", "--match-overlap", "--cis-normalization", "--trans-normalization", "--row-standardization", "--max-memory-gb", "--min-exposure-fraction", "-m", "--matrix-window-width", "-w", "--window", "--seed", "--compare", "--type", "--corr"};
                if (std::find(known.begin(), known.end(), name) == known.end()) throw std::invalid_argument("Unknown option: " + name);
                if (arg.find('=') != std::string::npos) value = arg.substr(arg.find('=') + 1);
                else if (++i < argc) value = argv[i];
                else throw std::invalid_argument("Missing value for " + name);
                if(name=="--k-search"){if(value!="screened"&&value!="exhaustive")throw std::invalid_argument("K search must be screened or exhaustive");options.kSearch=value;}
                else if(name=="--k-selection"){if(value!="bic"&&value!="aic")throw std::invalid_argument("K selection must be bic or aic");options.kSelection=value;options.adaptive=true;}
                else if(name=="--match-overlap"){size_t used=0;options.matchOverlap=std::stod(value,&used);if(used!=value.size()||!std::isfinite(options.matchOverlap)||options.matchOverlap<=.5||options.matchOverlap>1)throw std::invalid_argument("Match overlap must be in (0.5,1]");}
                else if(name=="--cis-view-weight"){size_t used=0;options.cisViewWeight=std::stod(value,&used);if(used!=value.size()||!std::isfinite(options.cisViewWeight)||options.cisViewWeight<0||options.cisViewWeight>100)throw std::invalid_argument("Cis view weight must be 0..100");}
                else if(name=="--normalization-views"){if(value!="inter"&&value!="dual")throw std::invalid_argument("views: inter or dual");options.dual=value=="dual";}
                else if(name=="--cluster-policy"){if(value!="adaptive"&&value!="fixed")throw std::invalid_argument("cluster policy: adaptive or fixed");options.adaptive=value=="adaptive";}
                else if(name=="--model"){if(value!="gmm"&&value!="kmeans")throw std::invalid_argument("adaptive model: gmm or kmeans");options.model=value;}
                else if(name=="--max-iterations")options.maxIterations=positive(value);
                else if(name=="--target-clusters")options.targetK=positive(value);
                else if(name=="--k-range"){auto range=split(value);if(range.size()!=2)throw std::invalid_argument("Expected MIN,MAX");options.adaptive=true;options.minK=positive(range[0]);options.maxK=positive(range[1]);if(options.maxK<options.minK||options.maxK>128)throw std::invalid_argument("Invalid K range");}
                else if(name=="--row-standardization") {if(value!="none"&&value!="weighted-zscore")throw std::invalid_argument("Row standardization must be none or weighted-zscore");options.rowStandardization=value!="none";}
                else if(name=="--max-memory-gb"||name=="--min-exposure-fraction") {size_t used=0;double v=std::stod(value,&used);if(used!=value.size()||!std::isfinite(v)||v<=0||(name=="--min-exposure-fraction"&&v>1))throw std::invalid_argument("Invalid exposure or memory limit");if(name=="--max-memory-gb")options.memoryGB=v;else options.minExposureFraction=v;}
                else if(name=="--trans-normalization") {options.transNormalization=value;options.balanceTrans=value=="NONE"||value=="AUTO";options.legacyCleaning=false;}
                else if(name=="--cis-normalization") options.normalization=value;
                else if (name == "--max-retries") options.maxRetries = positive(value);
                else if (name == "-r" || name == "--resolutions" || name == "--resolution") {
                    auto values = split(value); options.resolution = positive(values[0]);
                    if (values.size() > 1) std::cerr << "Only one resolution can be specified; using " << options.resolution << '\n';
                } else if (name == "-c" || name == "--chromosomes") options.chromosomes = split(value);
                else if (name == "-z" || name == "--threads") { auto n = integer(value); if (n > std::numeric_limits<int>::max()) throw std::invalid_argument("Thread count too large"); options.threads = n < 0 ? std::max(1U, std::thread::hardware_concurrency()) : std::max(1U, static_cast<unsigned>(n)); }
                else if (name == "-k" || name == "--normalization") {
                    if (value.empty()) throw std::invalid_argument("Normalization must not be empty");
                    options.normalization = value; options.transNormalization=value; options.balanceTrans=false; options.legacyCleaning=true;
                }
                else if (name == "--random-seeds") for (const auto& seed : split(value)) options.seed = static_cast<uint64_t>(integer(seed));
                else if (name == "--subsample") { auto n = integer(value); if (n > std::numeric_limits<int>::max() || n < std::numeric_limits<int>::min()) throw std::invalid_argument("Subsample out of 32-bit range"); if (n > 0) options.subsample = static_cast<int>(n); }
                else { std::cerr << "Warning: " << name << " is unused by Java slice and has no effect\n"; }
            } else positional.push_back(arg);
        }
        if (!positional.empty() && positional[0] == "slice") positional.erase(positional.begin());
        if (positional.empty()) { usage(); return argc == 1 ? 0 : 1; }
        if(positional.size()==2||positional.size()==3){
            options.files=split(positional[0]);options.output=positional[1];options.prefix=positional.size()==3?positional[2]:"slice";
            options.adaptive=true;options.firstK=1;options.endK=11;
            if(!options.minK){options.minK=1;options.maxK=10;}
        }else if(positional.size()==4){
            options.files=split(positional[0]);auto k=split(positional[1]);
            if(k.size()!=3)throw std::invalid_argument("Cluster argument must be K0,KF,nK");
            options.firstK=positive(k[0]);options.endK=positive(k[1]);options.attempts=positive(k[2]);
            if(options.endK<=options.firstK)throw std::invalid_argument("KF must be greater than K0");
            options.output=positional[2];options.prefix=positional[3];
        }else throw std::invalid_argument("Expected FILE.hic OUTPUT [PREFIX] or FILE.hic K0,KF,nK OUTPUT PREFIX");
        if (options.prefix.find('/') != std::string::npos || options.prefix.find('\\') != std::string::npos) throw std::invalid_argument("Prefix must be a filename prefix");
        if(options.forceBalance){options.balanceTrans=true;options.legacyCleaning=false;}
        if((options.dual||options.adaptive||options.model=="gmm")&&(options.javaCompatible||options.legacyCleaning))throw std::invalid_argument("Dual/adaptive fitting requires modern cleaning; omit legacy/Java switches");
        if(options.javaCompatible&&(options.rowStandardization||options.balanceTrans))throw std::invalid_argument("Row standardization is incompatible with --java-compatible");
        slice::run(options);
        return 0;
    } catch (const std::exception& e) { std::cerr << "slice: " << e.what() << '\n'; return 1; }
}
