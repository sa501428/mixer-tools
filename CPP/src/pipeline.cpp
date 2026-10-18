#include "slice/slice.hpp"
#include "pipeline_internal.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>

namespace slice {
namespace {
std::vector<int> labels(const std::vector<Cluster>& groups, size_t n) {
    std::vector<int> ids(n,-1);
    for (size_t c=0;c<groups.size();++c) for (auto i:groups[c].members) ids[i]=static_cast<int>(c);
    return ids;
}
double score(const Matrix& data, const std::vector<Cluster>& groups, bool medians) {
    double sum=0;
    for (const auto& c:groups) {
        if (c.members.size()<5) sum+=std::numeric_limits<float>::max();
        for (size_t i:c.members) sum+= medians?absoluteDistance(c.center,data[i]):squaredDistance(c.center,data[i]);
    }
    return sum/groups.size();
}
void plot(const std::filesystem::path& path, const Matrix& coordinates, const std::vector<int>& ids) {
    auto out=detail::outputFile(path);
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1000\" height=\"1000\" viewBox=\"0 0 1000 1000\"><rect width=\"1000\" height=\"1000\" fill=\"white\"/>\n";
    double lo[2]={1e100,1e100}, hi[2]={-1e100,-1e100};
    for (const auto& row:coordinates) for (int c=0;c<2;++c) { lo[c]=std::min(lo[c],double(row[c])); hi[c]=std::max(hi[c],double(row[c])); }
    for (size_t i=0;i<coordinates.size();++i) {
        double x=25+950*(coordinates[i][0]-lo[0])/std::max(1e-10,hi[0]-lo[0]);
        double y=25+950*(coordinates[i][1]-lo[1])/std::max(1e-10,hi[1]-lo[1]);
        auto c=ids[i]+1; // Java projection labels start at zero.
        // Stable categorical palette, independent of the embedding.
        out << "<circle cx=\""<<x<<"\" cy=\""<<y<<"\" r=\"2\" fill=\"hsl("<<(c*137)%360<<",65%,45%)\"/>\n";
    }
    out << "</svg>\n";
}
void encode(const Options& o, const detail::Input& input, const CleanMatrix& data,
            const std::vector<std::vector<int>>& all, const std::vector<std::vector<double>>& stats) {
    int best=-1;
    for (int z=int(all.size())-1;z>=0;--z) if (stats[5][z]>0 && stats[5][z]<.8) { best=z; break; }
    if (o.firstK<=5 && o.endK>5) best=std::max(best,5-o.firstK);
    if (best<0) throw std::runtime_error("ENCODE selection has no qualifying K (correlation in (0,.8) or K=5); expand the K range");
    auto out=detail::outputFile(o.output/"slice_subcompartment_clusters.bed");
    out << "#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRGB\tnumAltClusterings\taltClusterNum\taltClusterAssignment\n"
        << "#subcompartments called by mixer-tools/slice C++\n";
    for (size_t i=0;i<data.intervals.size();) {
        size_t end=i+1;
        while (end<data.intervals.size() && data.intervals[end].chromosome==data.intervals[i].chromosome && data.intervals[end].start==data.intervals[end-1].end) {
            bool same=true; for (const auto& ids:all) if (ids[end]!=ids[i]) same=false;
            if (!same) break;
            ++end;
        }
        // bedLine includes newline; write extended BED using a temporary stream.
        std::ostringstream line;
        detail::bedLine(line,input.chromosomes[data.intervals[i].chromosome],data.intervals[i].start,data.intervals[end-1].end,all[best][i]+1);
        auto text=line.str(); text.pop_back(); out<<text<<'\t'<<all.size()<<'\t';
        // Preserve the Java ENCODE writer's historical first-K+1 and trailing
        // width/zero fields, rather than silently changing its file schema.
        out<<o.firstK+1; for (size_t z=1;z<all.size();++z) out<<','<<o.firstK+z;
        out<<','<<data.intervals[end-1].end-data.intervals[i].start<<'\t';
        out<<all[0][i]; for (size_t z=1;z<all.size();++z) out<<','<<all[z][i];
        out<<",0\n"; i=end;
    }
}
}
void run(const Options& options) {
    auto inputs=detail::loadInputs(options);
    std::filesystem::create_directories(options.output);
    auto bad=detail::badBins(inputs,options.resolution);
    if (options.detectTranslocations) {
        for (auto& input:inputs) {
            int detection=0;
            input.translocations=detail::findTranslocations(input,bad,options.resolution,detection,options.output);
        }
    }
    Random command(options.seed), genome(command.nextLong());
    CleanMatrix data;
    uint64_t firstMatrixState=0;
    for (size_t d=0;d<inputs.size();++d) {
        Random composite(genome.nextLong());
        Random orderer(composite.nextLong());
        std::vector<detail::Order> orders;
        auto init=detail::outputFile(options.output/(d?"initial_split_"+std::to_string(d)+".bed":"initial_split.bed"));
        auto problems=detail::outputFile(options.output/(d?"problems_"+std::to_string(d)+".bed":"problems.bed"));
        for (size_t a=0;a<inputs[d].chromosomes.size();++a) {
            std::cout<<"Compressing "<<inputs[d].file<<" chromosome "<<inputs[d].chromosomes[a].name<<'\n';
            orders.push_back(detail::orderChromosome(inputs[d],a,bad[a],options,orderer));
            for (size_t x=0;x<orders.back().assignments.size();++x) {
                int id=orders.back().assignments[x];
                detail::bedLine(id<0?problems:init,inputs[d].chromosomes[a],int64_t(x)*options.resolution,int64_t(x+1)*options.resolution,id);
            }
        }
        auto block=detail::compressedTrans(inputs[d],bad,orders,options.resolution);
        if (!d) { data=std::move(block); composite.nextLong(); firstMatrixState=composite.nextLong(); }
        else {
            if (block.intervals.size()!=data.intervals.size()) throw std::runtime_error("Input row alignment mismatch");
            for (size_t i=0;i<data.values.size();++i) data.values[i].insert(data.values[i].end(),block.values[i].begin(),block.values[i].end());
            data.weights.insert(data.weights.end(),block.weights.begin(),block.weights.end());
        }
    }
    clean(data,options.filterOutliers);
    std::vector<std::vector<int32_t>> indices;
    for (const auto& interval:data.intervals) {
        if (interval.end>std::numeric_limits<int32_t>::max()) throw std::runtime_error("Java-compatible genome_indices.npy requires 32-bit coordinates");
        indices.push_back({inputs[0].chromosomes[interval.chromosome].index,static_cast<int32_t>(interval.start),static_cast<int32_t>(interval.end)});
    }
    detail::saveNpy(options.output/"genome_indices.npy",indices);
    detail::saveNpy(options.output/"slice_matrix.npy",data.values);
    Matrix coordinates;
    if (options.project) {
        auto features=cosineFeatures(data.values,firstMatrixState,options.threads);
        detail::saveNpy(options.output/"corr_slice_matrix.npy",features);
        std::cout<<"Generating diagnostic UMAP projection\n";
        coordinates=umap(features,options.threads);
        detail::saveNpy(options.output/"umap_projection.npy",coordinates);
        std::vector<int> chrom; for (auto interval:data.intervals) chrom.push_back(inputs[0].chromosomes[interval.chromosome].index);
        plot(options.output/"umap_chrom.svg",coordinates,chrom);
    }
    for (bool medians:{false,true}) {
        if (medians?!options.medians:!options.means) continue;
        if (medians) scaleWeights(data);
        std::string stem=medians?"kmedians":"kmeans";
        size_t nk=static_cast<size_t>(options.endK-options.firstK);
        std::vector<std::vector<double>> stats(6,std::vector<double>(nk));
        std::vector<std::vector<int>> all;
        for (int k=options.firstK;k<options.endK;++k) {
            double best=std::numeric_limits<double>::infinity(); std::vector<Cluster> winner;
            for (int attempt=0;attempt<options.attempts;++attempt) {
                std::vector<Cluster> groups;
                for (int retry=0;retry<options.maxRetries;++retry) {
                    groups=kmeans(data.values,k,200,genome.nextLong(),options.threads,medians);
                    if (groups.size()==static_cast<size_t>(k)) break;
                }
                if (groups.size()!=static_cast<size_t>(k)) throw std::runtime_error("Cannot obtain K="+std::to_string(k)+" nonempty clusters within retry limit");
                std::sort(groups.begin(),groups.end(),[](const Cluster& a,const Cluster& b) {
                    if (a.members.size()!=b.members.size()) return a.members.size()>b.members.size();
                    return a.members.front()<b.members.front();
                });
                double value=score(data.values,groups,medians);
                if (value<best) { best=value; winner=std::move(groups); }
            }
            if (winner.empty()) throw std::runtime_error("Clustering produced nonfinite scores");
            size_t z=static_cast<size_t>(k-options.firstK); double worst=-1;
            for (size_t a=0;a<winner.size();++a) for (size_t b=a+1;b<winner.size();++b) {
                double r=correlation(winner[a].center,winner[b].center);
                worst=std::isnan(r)?r:std::max(worst,r);
            }
            stats[0][z]=k; stats[1][z]=best; stats[2][z]=best+2*data.weights.size()*k;
            stats[3][z]=best+.5*data.weights.size()*k*std::log(data.values.size()); stats[4][z]=0; stats[5][z]=worst;
            auto ids=labels(winner,data.values.size()); all.push_back(ids);
            auto name=options.prefix+"_"+std::to_string(k)+"_"+stem+"_clusters";
            detail::saveBed(options.output/(name+".bed"),inputs[0],data,winner);
            if (!coordinates.empty()) plot(options.output/("umap_"+name+".svg"),coordinates,ids);
            std::cout<<stem<<" K="<<k<<" score="<<best<<'\n';
        }
        detail::saveNpy(options.output/(stem+"_cluster_post_stats.npy"),stats);
        if (options.encode && medians) encode(options,inputs[0],data,all,stats);
    }
    std::cout<<"Clustering complete\n";
}
}
