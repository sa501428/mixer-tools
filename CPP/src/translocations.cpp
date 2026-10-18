// Port of the dormant SimpleTranslocationFinder / Connected2DComponent helpers.
#include "pipeline_internal.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace slice::detail {
std::vector<Translocation> findTranslocations(const Input& input, const Mask& bad, int resolution,
    int& detectionResolution, const std::filesystem::path& output) {
    detectionResolution=*std::max_element(input.info.bpResolutions.begin(),input.info.bpResolutions.end());
    const int coarse=detectionResolution;
    const auto& norm=input.gw;
    for (const auto& chr:input.chromosomes) {
        std::vector<double> vector;
        if (!getNormalizationVectorForFile(input.file,chr.name,coarse,norm,vector) || vector.empty())
            throw std::runtime_error("Missing "+norm+" at the translocation detection resolution");
    }
    Mask masked;
    for (size_t a=0;a<bad.size();++a) {
        masked.emplace_back(bins(input.chromosomes[a],coarse),false);
        for (size_t x=0;x<bad[a].size();++x) if (bad[a][x]) masked[a][int64_t(x)*resolution/coarse]=true;
    }
    double mean=0,m2=0; size_t count=0;
    for (size_t a=0;a<bad.size();++a) for (size_t b=a+1;b<bad.size();++b)
        records(input,a,b,coarse,norm,[&](size_t x,size_t y,float v) {
            if (masked[a][x]||masked[b][y]||!std::isfinite(v)) return;
            double value=std::log1p(v),delta=value-mean;
            mean+=delta/++count; m2+=delta*(value-mean);
        });
    auto suffix=std::to_string(coarse)+"_"+std::filesystem::path(input.file).filename().string()+".bedpe";
    auto all=outputFile(output/("all_suspected_translocations_"+suffix));
    auto selected=outputFile(output/("suspected_translocations_"+suffix));
    std::vector<Translocation> result;
    if (count<2||!(m2>0)) return result;
    double sd=std::sqrt(m2/(count-1));
    for (size_t a=0;a<bad.size();++a) for (size_t b=a+1;b<bad.size();++b) {
        struct Point {size_t x,y;double z;}; std::vector<Point> points;
        size_t rows=0,cols=0;
        records(input,a,b,coarse,norm,[&](size_t x,size_t y,float v) {
            if (masked[a][x]||masked[b][y]||!std::isfinite(v)) return;
            double z=(std::log1p(v)-mean)/sd;
            if (z>6) {
                points.push_back({x,y,z}); rows=std::max(rows,x+1); cols=std::max(cols,y+1);
                all<<input.chromosomes[a].name<<'\t'<<x*coarse<<'\t'<<(x+1)*coarse<<'\t'<<input.chromosomes[b].name<<'\t'<<y*coarse<<'\t'<<(y+1)*coarse<<'\t'<<z<<'\n';
            }
        });
        if (points.empty()) continue;
        if (rows>100000000/cols) throw std::runtime_error("Translocation bitmap exceeds 100 million cells; use a coarser map");
        Mask bitmap(rows,std::vector<bool>(cols,false));
        for (auto p:points) for (size_t x=p.x?p.x-1:0;x<std::min(rows,p.x+2);++x)
            for (size_t y=p.y?p.y-1:0;y<std::min(cols,p.y+2);++y) bitmap[x][y]=true;
        for (size_t x=0;x<rows;++x) for (size_t y=0;y+1<cols;++y) if (bitmap[x][y]&&!bitmap[x][y+1])
            for (size_t end=std::min(cols-1,y+10);end>y;--end) if (bitmap[x][end]) {
                for (size_t j=y+1;j<end;++j) bitmap[x][j]=true;
                break;
            }
        for (size_t x=0;x+1<rows;++x) for (size_t y=0;y<cols;++y) if (bitmap[x][y]&&!bitmap[x+1][y])
            for (size_t end=std::min(rows-1,x+10);end>x;--end) if (bitmap[end][y]) {
                for (size_t j=x+1;j<end;++j) bitmap[j][y]=true;
                break;
            }
        // Preserve Java's bounding-box expansion, including filling enclosed
        // disconnected islands, rather than substituting flood-fill components.
        for (size_t x=0;x<rows;++x) for (size_t y=0;y<cols;++y) if (bitmap[x][y]) {
            int x0=int(x),x1=x0+1,y0=int(y),y1=y0+1;
            auto any=[&](int u0,int u1,int v0,int v1) {
                for (int u=std::max(0,u0);u<std::min(int(rows),u1);++u)
                    for (int v=std::max(0,v0);v<std::min(int(cols),v1);++v) if (bitmap[u][v]) return true;
                return false;
            };
            bool change=true;
            while (change) {
                change=false;
                if (any(x0,x1+1,y1,y1+1)) {++y1;change=true;}
                if (any(x0,x1+1,y0-1,y0)) {--y0;change=true;}
                if (any(x1,x1+1,y0-1,y1+1)) {++x1;change=true;}
            }
            for (int u=x0;u<x1;++u) for (int v=y0;v<y1;++v) bitmap[u][v]=false;
            if ((x1-x0)*(y1-y0)<=20) continue;
            double peak=0; for (auto p:points) if (int(p.x)>=x0&&int(p.x)<x1&&int(p.y)>=y0&&int(p.y)<y1) peak=std::max(peak,p.z);
            result.push_back({a,b,int64_t(x0)*coarse,int64_t(x1)*coarse,int64_t(y0)*coarse,int64_t(y1)*coarse,peak});
            selected<<input.chromosomes[a].name<<'\t'<<int64_t(x0)*coarse<<'\t'<<int64_t(x1)*coarse<<'\t'<<input.chromosomes[b].name<<'\t'<<int64_t(y0)*coarse<<'\t'<<int64_t(y1)*coarse<<'\n';
        }
    }
    return result;
}
}
