#pragma once
// Shared BED9 presentation contract. Keep copies in the caller repositories identical.
#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace subcompartment_bed {
inline std::string color(int64_t id) {
    // High-contrast categorical colors (Kelly-inspired); no light/dark pairs.
    static const char* palette[] = {
        "0,114,178", "230,159,0", "0,158,115", "213,62,0", "170,0,170",
        "0,190,210", "240,80,150", "110,70,30", "190,190,0", "50,40,120",
        "140,210,40", "35,35,35", "255,180,160", "100,0,45", "0,90,70",
        "180,140,220", "240,210,50", "100,160,230", "150,100,100", "0,45,65"
    };
    if (id <= 0) return "128,128,128"; // Explicit unassigned state.
    if (id <= 20) return palette[id-1];
    // An odd multiplier permutes 24-bit RGB values instead of cycling at K=20.
    uint32_t rgb = (uint64_t(id-20)*0x9e3779u) & 0xffffffu;
    return std::to_string((rgb>>16)&255)+","+std::to_string((rgb>>8)&255)+","+std::to_string(rgb&255);
}
struct Interval { std::string chrom, label; int64_t start, end; };
inline void line(std::ostream& out, const Interval& x, const std::string& rgb) {
    out << x.chrom << '\t' << x.start << '\t' << x.end << '\t' << x.label
        << "\t0\t.\t" << x.start << '\t' << x.end << '\t' << rgb << '\n';
}
// Explicitly invoked only for call tracks, after their producer has closed them.
// Preserve chromosome names, labels and coverage; only sort within chromosomes,
// merge exactly touching equal labels and replace display columns with BED9.
inline void finish(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot read call BED: " + path.string());
    std::vector<Interval> rows;
    std::map<std::string,size_t> chromOrder;
    std::set<std::string> labels;
    std::string text;
    while (std::getline(in,text)) {
        if (text.empty() || text[0]=='#' || text.rfind("track ",0)==0) continue;
        std::istringstream stream(text); Interval x;
        if (!(stream >> x.chrom >> x.start >> x.end >> x.label) || x.start<0 || x.end<=x.start)
            throw std::runtime_error("Invalid call BED row: " + path.string());
        if (!chromOrder.count(x.chrom)) chromOrder[x.chrom]=chromOrder.size();
        rows.push_back(x); labels.insert(x.label);
    }
    if (in.bad()) throw std::runtime_error("Cannot read call BED: " + path.string());
    in.close();
    std::stable_sort(rows.begin(),rows.end(),[&](const Interval& a,const Interval& b){
        return std::make_tuple(chromOrder.at(a.chrom),a.start,a.end)<std::make_tuple(chromOrder.at(b.chrom),b.start,b.end);
    });
    std::map<std::string,std::string> colors;
    std::set<std::string> used;
    int64_t next=1;
    for (const auto& label:labels) {
        size_t n=0; int64_t id=0;
        try {id=std::stoll(label,&n);} catch (const std::exception&) {}
        if(n==label.size()) {colors[label]=color(id);used.insert(colors[label]);}
    }
    for (const auto& label:labels) if(!colors.count(label)) {
        std::string rgb; do {rgb=color(next++);} while(used.count(rgb));
        colors[label]=rgb;used.insert(rgb);
    }
    auto temp=path;temp += ".bed9.tmp";
    std::ofstream out(temp);
    if(!out) throw std::runtime_error("Cannot write call BED: "+temp.string());
    for(size_t i=0;i<rows.size();) {
        Interval x=rows[i++];
        while(i<rows.size() && rows[i].chrom==x.chrom && rows[i].start==x.end && rows[i].label==x.label)
            x.end=rows[i++].end;
        line(out,x,colors.at(x.label));
    }
    out.close();
    if(!out) throw std::runtime_error("Cannot finish call BED: "+temp.string());
    std::filesystem::rename(temp,path);
}
} // namespace subcompartment_bed
