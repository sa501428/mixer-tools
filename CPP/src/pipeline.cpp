#include "slice/bed_output.hpp"
#include "slice/slice.hpp"
#include "slice/trans_preprocess.hpp"
#include "slice/adaptive.hpp"
#include "slice/parallel.hpp"
#include <chrono>
#include <map>
#include <mutex>
#include <set>
#include <tuple>
#include "pipeline_internal.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <numeric>

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
    auto out=detail::outputFile(o.output/"slice_subcompartment_clusterings.tsv");
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
    auto bed=detail::outputFile(o.output/"slice_subcompartment_clusters.bed");
    for(size_t i=0;i<data.intervals.size();++i) {
        const auto& x=data.intervals[i];
        detail::bedLine(bed,input.chromosomes[x.chromosome],x.start,x.end,all[best][i]+1);
    }
    bed.close();
    subcompartment_bed::finish(o.output/"slice_subcompartment_clusters.bed");
}
}
namespace {
void adaptiveRun(const Options& o,const detail::Input& input,const CleanMatrix& inter,const CleanMatrix* cis){
    int target=o.targetK?o.targetK:o.firstK;if(target<1||target>128)throw std::runtime_error("Adaptive target K must be 1..128");int lo=o.minK?o.minK:(o.adaptive?std::max(1,target-2):target),hi=o.maxK?o.maxK:(o.adaptive?target+4:target);
    if(target<1||target>128||hi>128)throw std::runtime_error("Adaptive maximum K exceeds 128");
    std::vector<const CleanMatrix*> views{&inter};if(cis)views.push_back(cis);
    std::vector<hic_adaptive::Fit> fits;auto stats=detail::outputFile(o.output/"adaptive_fits.tsv");stats<<"view\ttarget_k\tselected_k\tbic\tconverged\titerations\n";
    for(size_t v=0;v<views.size();++v){const auto& data=*views[v];hic_adaptive::Rows rows;for(const auto& row:data.values)rows.emplace_back(row.begin(),row.end());
        struct Job{int k,attempt;};std::vector<Job> jobs;for(int k=lo;k<=hi&&size_t(k*3)<=rows.size();++k)for(int a=0;a<o.attempts;++a)jobs.push_back({k,a});
        // Independent starts share read-only features, with bounded posterior workspaces.
        double bytes=std::max(1.,double(rows.size())*hi*16);unsigned workers=std::max(1u,std::min(o.threads,unsigned(std::max(1.,o.memoryGB*1073741824*.2/bytes))));
        workers=std::max(1u,std::min(workers,unsigned(jobs.size())));
        std::mutex resultMutex;hic_adaptive::Fit best;size_t bestIndex=jobs.size();
        detail::parallelFor(jobs.size(),workers,[&](size_t j){const auto& job=jobs[j];auto fit=hic_adaptive::fit(rows,job.k,o.maxIterations,o.seed+1009*job.k+104729*job.attempt,o.model=="gmm",{},std::max(1u,o.threads/workers));
            std::lock_guard<std::mutex> lock(resultMutex);if(fit.bic<best.bic||(fit.bic==best.bic&&j<bestIndex)){best=std::move(fit);bestIndex=j;}});
        if(!std::isfinite(best.bic)){if(v){std::cerr<<"Warning: no supported cis fit; returning inter calls\n";views.resize(1);break;}throw std::runtime_error("No supported adaptive Mixer fit");}
        stats<<(v?"cis":"inter")<<'\t'<<target<<'\t'<<best.means.size()<<'\t'<<best.bic<<'\t'<<best.converged<<'\t'<<best.iterations<<'\n';
        std::vector<Cluster> clusters(best.means.size());for(size_t i=0;i<best.labels.size();++i)if(best.labels[i]>=0)clusters[best.labels[i]].members.push_back(i);
        detail::saveBed(o.output/(v?"cis.bed":"inter.bed"),input,data,clusters);
        auto probabilities=detail::outputFile(o.output/(v?"cis_memberships.tsv":"inter_memberships.tsv"));probabilities<<"chromosome\tstart\tend\tlabel\tsoft_score\n";
        for(size_t i=0;i<data.intervals.size();++i){auto x=data.intervals[i];probabilities<<input.chromosomes[x.chromosome].name<<'\t'<<x.start<<'\t'<<x.end<<'\t'<<best.labels[i]+1<<'\t'<<*std::max_element(best.probability[i].begin(),best.probability[i].end())<<'\n';}
        auto model=detail::outputFile(o.output/(v?"cis_model.tsv":"inter_model.tsv"));model<<"component\tfeature\tmean\tvariance\tprior\n";for(size_t c=0;c<best.means.size();++c)for(size_t j=0;j<best.means[c].size();++j)model<<c+1<<'\t'<<j+1<<'\t'<<best.means[c][j]<<'\t'<<best.variance[c][j]<<'\t'<<best.prior[c]<<'\n';
        fits.push_back(std::move(best));
    }
    using Key=std::tuple<size_t,int64_t,int64_t>;std::map<Key,size_t> index;std::vector<Interval> intervals;
    for(auto view:views)for(const auto& x:view->intervals){Key key{x.chromosome,x.start,x.end};if(!index.count(key)){index[key]=intervals.size();intervals.push_back(x);}}
    std::vector<std::vector<int>> tracks(views.size(),std::vector<int>(index.size(),-1));for(size_t v=0;v<views.size();++v)for(size_t i=0;i<views[v]->intervals.size();++i){auto x=views[v]->intervals[i];tracks[v][index.at({x.chromosome,x.start,x.end})]=fits[v].labels[i];}
    auto mapping=hic_adaptive::reconcile(tracks,{},o.matchOverlap);auto calls=hic_adaptive::consensus(tracks,mapping,tracks.size()>1?std::vector<double>{1,o.cisViewWeight}:std::vector<double>{1});hic_adaptive::lineage((o.output/"state_lineage.tsv").string(),mapping);
    auto bed=detail::outputFile(o.output/"combined.bed");size_t called=0,both=0;std::vector<std::pair<Key,size_t>> sorted(index.begin(),index.end());std::set<int> states;
    auto uncertainty=detail::outputFile(o.output/"ambiguous_bins.tsv");uncertainty<<"chromosome\tstart\tend\tinter_state\tcis_state\n";
    auto footprint=detail::outputFile(o.output/"view_footprint.tsv");footprint<<"chromosome\tstart\tend\tinter_available\tcis_available\n";
    for(const auto& item:sorted){size_t i=item.second;auto x=intervals[i];if(tracks.size()>1&&hic_adaptive::unresolved(tracks[0][i],tracks[1][i],mapping))uncertainty<<input.chromosomes[x.chromosome].name<<'\t'<<x.start<<'\t'<<x.end<<'\t'<<tracks[0][i]+1<<'\t'<<(tracks.size()>1?tracks[1][i]+1:0)<<'\n';if(calls[i]>=0){detail::bedLine(bed,input.chromosomes[x.chromosome],x.start,x.end,calls[i]+1);states.insert(calls[i]);called++;}bool c=tracks.size()>1&&tracks[1][i]>=0;both+=tracks[0][i]>=0&&c;footprint<<input.chromosomes[x.chromosome].name<<'\t'<<x.start<<'\t'<<x.end<<'\t'<<(tracks[0][i]>=0)<<'\t'<<c<<'\n';}
    bed.close();
    subcompartment_bed::finish(o.output/"combined.bed");
    auto report=detail::outputFile(o.output/"adaptive_summary.tsv");report<<"target_k\treturned_k\tcandidate_states\tcalled_bins\tshared_bins\n"<<target<<'\t'<<states.size()<<'\t'<<mapping.members.size()<<'\t'<<called<<'\t'<<both<<'\n';
    std::cout<<"Target "<<target<<"; returned "<<states.size()<<" states; "<<std::count(mapping.provisional.begin(),mapping.provisional.end(),true)<<" provisional candidates\n";
}
}
void run(const Options& options) {
    auto inputs=detail::loadInputs(options);
    std::filesystem::create_directories(options.output);
    auto parameters=detail::outputFile(options.output/"analysis_parameters.tsv");
    parameters<<"parameter\tvalue\ntrans_mode\t"<<(options.balanceTrans?"computed_if_raw":"stored")<<"\nrow_standardization\t"<<(options.rowStandardization?"weighted-zscore":"none")<<"\ncleaning\t"<<(options.legacyCleaning?"legacy":"zero_preserving_log1p_clip")<<"\nmax_memory_gb\t"<<options.memoryGB<<'\n';
    parameters<<"cluster_policy\t"<<(options.adaptive?"adaptive":"fixed")<<"\nnormalization_views\t"<<(options.dual?"dual":"inter")<<"\ncis_view_weight\t"<<options.cisViewWeight<<"\nmatch_overlap\t"<<options.matchOverlap<<"\nmodel\t"<<((options.adaptive||options.dual||options.model=="gmm")?options.model:(options.means&&options.medians?"kmeans,kmedians":options.medians?"kmedians":"kmeans"))<<"\nmax_iterations\t"<<options.maxIterations<<"\nthreads\t"<<options.threads<<"\ndistance_bounds\t"<<options.distanceBounds<<'\n';
    parameters<<"cis_coverage_statistic\traw_off_diagonal_nonzeros\ncis_coverage_z\t-2.5\ncis_coverage_population\tnonempty_rows_per_chromosome\n";
    auto bad=detail::badBins(inputs,options.resolution,{},options.output/"cis_coverage.tsv");
    if (options.detectTranslocations) {
        for (auto& input:inputs) {
            int detection=0;
            input.translocations=detail::findTranslocations(input,bad,options.resolution,detection,options.output,options.memoryGB);
        }
    }
    Random command(options.seed), genome(command.nextLong());
    CleanMatrix data,cisData;
    bool cisAvailable=options.dual;
    auto began=std::chrono::steady_clock::now();
    uint64_t firstMatrixState=0;
    for (size_t d=0;d<inputs.size();++d) {
        Random composite(genome.nextLong());
        Random orderer(composite.nextLong());
        std::vector<detail::Order> orders;
        auto init=detail::outputFile(options.output/(d?"initial_split_"+std::to_string(d)+".bed":"initial_split.bed"));
        auto problems=detail::outputFile(options.output/(d?"problems_"+std::to_string(d)+".bed":"problems.bed"));
        if(options.javaCompatible){for(size_t a=0;a<inputs[d].chromosomes.size();++a)orders.push_back(detail::orderChromosome(inputs[d],a,bad[a],options,orderer));}
        else {
        orders.resize(inputs[d].chromosomes.size());
        std::vector<Random> chromosomeRandom;double largest=1;
        for(const auto& chr:inputs[d].chromosomes){chromosomeRandom.push_back(orderer);orderer.nextLong();double n=double(chr.length/std::max(100000,options.resolution)+1);largest=std::max(largest,96*n*n);}
        unsigned workers=options.javaCompatible?1:std::max(1u,std::min(options.threads,unsigned(std::max(1.,options.memoryGB*1073741824*.6/largest))));
        workers=std::min<unsigned>(workers,unsigned(orders.size()));
        Options compression=options;compression.threads=std::max(1u,options.threads/workers);compression.memoryGB=options.memoryGB/workers;
        std::cout<<"Compressing "<<orders.size()<<" chromosomes with "<<workers<<" independent workers\n";
        detail::parallelFor(orders.size(),workers,[&](size_t a){orders[a]=detail::orderChromosome(inputs[d],a,bad[a],compression,chromosomeRandom[a]);});
        }
        for (size_t a=0;a<inputs[d].chromosomes.size();++a) {
            std::cout<<"Compressed "<<inputs[d].file<<" chromosome "<<inputs[d].chromosomes[a].name<<'\n';
            for (size_t x=0;x<orders[a].assignments.size();++x) {
                int id=orders[a].assignments[x];
                detail::bedLine(id<0?problems:init,inputs[d].chromosomes[a],int64_t(x)*options.resolution,int64_t(x+1)*options.resolution,id);
            }
        }
        auto block=detail::compressedTrans(inputs[d],bad,orders,options.resolution,options.minExposureFraction,options.memoryGB,options.legacyCleaning,options.threads);
        if(cisAvailable){auto cisInput=inputs[d];cisInput.inter=options.normalization.empty()?detail::chooseNormLabel(cisInput,{"SCALE","KR","NONE"}):options.normalization;
            if(cisInput.inter=="NONE"){std::cerr<<"Warning: cis companion unavailable with NONE; inter analysis continues\n";cisAvailable=false;}
            else {auto cb=detail::compressedTrans(cisInput,bad,orders,options.resolution,options.minExposureFraction,options.memoryGB,false,options.threads,true);
                if(!d)cisData=std::move(cb);else {hic_trans::memoryGuard(cb.values.size(),cisData.weights.size()+cb.weights.size(),options.memoryGB);for(size_t i=0;i<cb.values.size();++i){cisData.values[i].insert(cisData.values[i].end(),cb.values[i].begin(),cb.values[i].end());cisData.exposure[i].insert(cisData.exposure[i].end(),cb.exposure[i].begin(),cb.exposure[i].end());}cisData.weights.insert(cisData.weights.end(),cb.weights.begin(),cb.weights.end());}}
        }
        parameters<<"cis_normalization_"<<d<<'\t'<<inputs[d].intra<<"\ntrans_source_normalization_"<<d<<'\t'<<inputs[d].inter<<'\n';
        if(options.balanceTrans&&(inputs[d].inter=="NONE"||options.forceBalance)){
            auto report=hic_trans::balance(block);
            parameters<<"balance_iterations_"<<d<<'\t'<<report.iterations<<"\nbalance_error_"<<d<<'\t'<<report.relativeError<<"\nbalance_mean_contact_"<<d<<'\t'<<report.meanContact<<'\n';
            auto rf=detail::outputFile(options.output/("trans_row_factors_"+std::to_string(d)+".tsv"));rf<<"row\tchromosome\tstart\tend\tmultiplier\tavailable_target_mass\n";for(size_t i=0;i<report.rowFactors.size();++i){const auto& x=block.intervals[i];rf<<i<<'\t'<<inputs[d].chromosomes[x.chromosome].name<<'\t'<<x.start<<'\t'<<x.end<<'\t'<<report.rowFactors[i]<<'\t'<<std::accumulate(block.exposure[i].begin(),block.exposure[i].end(),0.)<<'\n';}
            auto cf=detail::outputFile(options.output/("trans_column_factors_"+std::to_string(d)+".tsv"));cf<<"column\tmultiplier\n";for(size_t j=0;j<report.columnFactors.size();++j)cf<<j<<'\t'<<report.columnFactors[j]<<'\n';
        }
        if (!d) { data=std::move(block); composite.nextLong(); firstMatrixState=composite.nextLong(); }
        else {
            if (block.intervals.size()!=data.intervals.size()) throw std::runtime_error("Input row alignment mismatch");
            hic_trans::memoryGuard(data.values.size(),data.weights.size()+block.weights.size(),options.memoryGB);
            for (size_t i=0;i<data.values.size();++i){ data.values[i].insert(data.values[i].end(),block.values[i].begin(),block.values[i].end());data.exposure[i].insert(data.exposure[i].end(),block.exposure[i].begin(),block.exposure[i].end());}
            data.weights.insert(data.weights.end(),block.weights.begin(),block.weights.end());
        }
    }
    if(options.legacyCleaning)clean(data,options.filterOutliers);else hic_trans::clean(data,options.rowStandardization,options.filterOutliers);
    std::cout<<"Compression/cleaning: "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<" seconds\n";
    if(options.dual||options.adaptive||options.model=="gmm"){
        if(cisAvailable){try{hic_trans::clean(cisData,options.rowStandardization,options.filterOutliers);}catch(const std::runtime_error& e){if(std::string(e.what())!="No rows have variable features after normalization")throw;std::cerr<<"Warning: no variable cis companion features; inter analysis continues\n";cisAvailable=false;}}
        adaptiveRun(options,inputs[0],data,cisAvailable?&cisData:nullptr);return;
    }
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
            bool parallelDone=false;
            double scratch=double(data.values.size())*(16+8*k)+(medians?4.*data.values.size()*data.weights.size():0);
            unsigned fitWorkers=std::max(1u,std::min({options.threads,unsigned(options.attempts),unsigned(std::max(1.,options.memoryGB*1073741824*.25/std::max(1.,scratch)))}));
            if(options.parallelAttempts&&fitWorkers>1){
                Random before=genome;std::vector<uint64_t> seeds;for(int a=0;a<options.attempts;++a)seeds.push_back(genome.nextLong());
                struct Attempt{std::vector<Cluster> groups;double value=INFINITY;};std::vector<Attempt> attempts(size_t(options.attempts));
                detail::parallelFor(attempts.size(),fitWorkers,[&](size_t a){unsigned workers=std::max(1u,options.threads/fitWorkers);
                    attempts[a].groups=kmeans(data.values,k,options.maxIterations,seeds[a],workers,medians,options.distanceBounds);
                    if(attempts[a].groups.size()==size_t(k)){auto& groups=attempts[a].groups;std::sort(groups.begin(),groups.end(),[](const Cluster&a,const Cluster&b){if(a.members.size()!=b.members.size())return a.members.size()>b.members.size();return a.members.front()<b.members.front();});attempts[a].value=score(data.values,groups,medians);}});
                parallelDone=std::all_of(attempts.begin(),attempts.end(),[&](const Attempt&a){return a.groups.size()==size_t(k);});
                if(parallelDone){for(auto& attempt:attempts)if(attempt.value<best){best=attempt.value;winner=std::move(attempt.groups);}}
                else genome=before; // Replay retries serially to preserve the original Java seed stream.
            }
            if(!parallelDone)
            for (int attempt=0;attempt<options.attempts;++attempt) {
                std::vector<Cluster> groups;
                for (int retry=0;retry<options.maxRetries;++retry) {
                    groups=kmeans(data.values,k,options.maxIterations,genome.nextLong(),options.threads,medians,options.distanceBounds);
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
