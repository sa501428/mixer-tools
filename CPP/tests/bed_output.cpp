#include "slice/bed_output.hpp"
#include <iostream>

int main(int argc, char** argv) {
    if(argc!=2) return 2;
    const std::filesystem::path path=argv[1];
    auto require=[](bool condition){if(!condition) throw std::runtime_error("BED output contract failed");};
    std::set<std::string> palette;
    for(int id=0;id<=128;++id) require(palette.insert(subcompartment_bed::color(id)).second);
    {
        std::ofstream out(path);
        // Out-of-order members, a gap, a label transition, chromosome boundary,
        // unassigned bins, string labels and a clipped final interval.
        out << "chr1\t10\t20\t1\nchr1\t0\t10\t1\nchr1\t30\t40\t1\n"
            << "chr1\t40\t50\t2\nchr1\t50\t53\t2\nchr2\t0\t10\t2\n"
            << "chr2\t10\t20\t0\nchr2\t20\t30\t0\nchr2\t30\t40\tA.1\nchr2\t40\t43\tA.1\n";
    }
    subcompartment_bed::finish(path);
    auto read=[&](){std::ifstream in(path);return std::string(std::istreambuf_iterator<char>(in),{});};
    const auto result=read();
    require(std::count(result.begin(),result.end(),'\n')==6);
    require(result.find("chr1\t0\t20\t1\t0\t.\t0\t20\t0,114,178\n")!=std::string::npos);
    require(result.find("chr1\t40\t53\t2\t0\t.\t40\t53\t230,159,0\n")!=std::string::npos);
    require(result.find("chr2\t10\t30\t0\t0\t.\t10\t30\t128,128,128\n")!=std::string::npos);
    require(result.find("chr2\t30\t43\tA.1\t0\t.\t30\t43\t")!=std::string::npos);
    subcompartment_bed::finish(path); require(read()==result);
    {std::ofstream out(path);}
    subcompartment_bed::finish(path);require(read().empty());
    std::filesystem::remove(path);
    std::cout<<"BED9, adjacency, gaps, chromosome boundaries, clipped ends, palette and idempotence passed\n";
}
