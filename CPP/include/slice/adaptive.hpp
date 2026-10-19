#pragma once
#include <algorithm>
#include <cmath>
#include <fstream>
#include <atomic>
#include <future>
#include <mutex>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// Shared by the three native callers. Missing coordinates never become zeros.
namespace hic_adaptive {
using Rows=std::vector<std::vector<double>>;
struct Fit {
    Rows means,variance,logVariance,probability;
    std::vector<double> prior;
    std::vector<int> labels;
    double bic=INFINITY,aic=INFINITY,likelihood=-INFINITY;
    double score(const std::string& criterion) const { return criterion=="aic"?aic:bic; }
    bool converged=false;
    int iterations=0;
};
inline void posteriorInto(const std::vector<double>& x,const Fit& f,bool gaussian,std::vector<double>& score) {
    score.resize(f.means.size());
    for(size_t c=0;c<score.size();++c){double v=0;size_t n=0;
        for(size_t j=0;j<x.size();++j)if(std::isfinite(x[j])&&std::isfinite(f.means[c][j])){
            double delta=x[j]-f.means[c][j],var=gaussian?f.variance[c][j]:1.;
            v+=delta*delta/var+(gaussian?(f.logVariance.empty()?std::log(6.283185307179586*var):f.logVariance[c][j]):0);++n;}
        score[c]=n?(gaussian?std::log(std::max(1e-12,f.prior[c])):0)-.5*v:-INFINITY;
    }
}
inline std::vector<double> posterior(const std::vector<double>& x,const Fit& f,bool gaussian){std::vector<double> result;posteriorInto(x,f,gaussian,result);return result;}
template<class Function> inline void parallelBlocks(size_t n,unsigned threads,Function function){
    if(threads<=1||n<1024){function(0,n);return;}
    std::atomic<size_t> next{0};auto work=[&]{for(size_t begin;(begin=next.fetch_add(128))<n;)function(begin,std::min(n,begin+128));};
    std::vector<std::future<void>> workers;for(size_t t=1;t<std::min<size_t>(threads,(n+127)/128);++t)workers.push_back(std::async(std::launch::async,work));work();for(auto& worker:workers)worker.get();
}
inline Fit fit(const Rows& x,int k,int iterations,uint64_t seed,bool gaussian=false,
               const std::vector<double>& observationWeight={},unsigned threads=1) {
    Fit f;if(x.empty()||k<1||size_t(k)>x.size())return f;
    size_t n=x.size(),d=x[0].size();if(!d)return f;
    if(std::any_of(x.begin(),x.end(),[](const auto& row){return std::none_of(row.begin(),row.end(),[](double v){return std::isfinite(v);});}))return f;
    if(!observationWeight.empty()&&observationWeight.size()!=n)throw std::invalid_argument("Observation weight dimensions");
    std::vector<double> weights=observationWeight.empty()?std::vector<double>(n,1):observationWeight;
    double neff=0,weight2=0;for(double w:weights){if(!(w>0&&std::isfinite(w)))throw std::invalid_argument("Invalid observation weight");neff+=w;weight2+=w*w;}
    // Normalize fitting mass to the number of independent observations.
    for(double& w:weights)w*=n/neff;neff=neff*neff/weight2;
    std::vector<double> global(d),squares(d),mass(d),var(d,1);
    for(const auto& row:x){if(row.size()!=d)throw std::invalid_argument("Ragged features");for(size_t j=0;j<d;++j)if(std::isfinite(row[j])){global[j]+=row[j];squares[j]+=row[j]*row[j];mass[j]++;}}
    size_t modeledFeatures=size_t(std::count_if(mass.begin(),mass.end(),[](double v){return v>0;}));
    for(size_t j=0;j<d;++j)global[j]=mass[j]?global[j]/mass[j]:NAN;
    for(size_t j=0;j<d;++j)var[j]=mass[j]?std::max(1e-10,squares[j]/mass[j]-global[j]*global[j]):1;
    std::mt19937_64 rng(seed);std::vector<double> nearest(n,1);
    for(int c=0;c<k;++c){std::discrete_distribution<size_t> pick(nearest.begin(),nearest.end());size_t r=pick(rng);auto center=x[r];
        // Initialization only: global means for unobserved center coordinates.
        for(size_t j=0;j<d;++j)if(!std::isfinite(center[j]))center[j]=global[j];
        f.means.push_back(center);f.variance.push_back(var);
        for(size_t i=0;i<n;++i){double s=0;size_t m=0;for(size_t j=0;j<d;++j)if(std::isfinite(x[i][j])&&std::isfinite(center[j])){s+=(x[i][j]-center[j])*(x[i][j]-center[j]);m++;}
            nearest[i]=std::min(c?nearest[i]:INFINITY,m?s/m:0.);}}
    f.prior.assign(k,1./k);f.probability.assign(n,std::vector<double>(k));f.labels.resize(n);
    double previous=-INFINITY;int repairs=0;
    for(int it=0;it<iterations;++it){double ll=0;std::vector<double> cm(k),loss(n);
        parallelBlocks(n,threads,[&](size_t begin,size_t end){std::vector<double> score;
            for(size_t i=begin;i<end;++i){posteriorInto(x[i],f,gaussian,score);double best=*std::max_element(score.begin(),score.end());if(!std::isfinite(best))continue;
                double sum=0;for(double value:score)sum+=std::exp(value-best);loss[i]=best+std::log(sum);
                size_t label=size_t(std::max_element(score.begin(),score.end())-score.begin());f.labels[i]=int(label);
                for(int c=0;c<k;++c)f.probability[i][c]=gaussian?std::exp(score[c]-best)/sum:double(size_t(c)==label);}});
        // Fixed row-order reductions retain deterministic models across thread counts.
        for(size_t i=0;i<n;++i){ll+=weights[i]*loss[i];for(int c=0;c<k;++c)cm[c]+=weights[i]*f.probability[i][c];}
        bool repaired=false;
        Rows sx(k,std::vector<double>(d)),ss=sx,featureMass=sx;
        for(size_t i=0;i<n;++i)for(int c=0;c<k;++c){double w=weights[i]*f.probability[i][c];if(w<1e-15)continue;
            for(size_t j=0;j<d;++j)if(std::isfinite(x[i][j])){featureMass[c][j]+=w;sx[c][j]+=w*x[i][j];if(gaussian)ss[c][j]+=w*x[i][j]*x[i][j];}}

        for(int c=0;c<k;++c){if(cm[c]<.5){if(repairs>=2*k)return Fit{};size_t worst=size_t(std::min_element(loss.begin(),loss.end())-loss.begin());f.means[c]=x[worst];for(size_t j=0;j<d;++j)if(!std::isfinite(f.means[c][j]))f.means[c][j]=global[j];f.variance[c]=var;f.prior[c]=1./k;loss[worst]=INFINITY;repairs++;repaired=true;continue;}
            f.prior[c]=cm[c]/n;
            for(size_t j=0;j<d;++j){double m=featureMass[c][j];f.means[c][j]=m>1e-8?sx[c][j]/m:global[j];
                if(gaussian)f.variance[c][j]=m>1e-8?std::max(1e-10,std::max(.01*var[j],ss[c][j]/m-f.means[c][j]*f.means[c][j]+.001*var[j])):var[j];}}
        if(gaussian){f.logVariance=f.variance;for(auto& row:f.logVariance)for(double& v:row)v=std::log(6.283185307179586*v);}
        f.iterations=it+1;
        if(!repaired&&it>3&&std::abs(ll-previous)<1e-6*(1+std::abs(previous))){f.converged=true;break;}previous=repaired?-INFINITY:ll;
    }
    f.likelihood=0;std::vector<size_t> support(k);
    for(size_t i=0;i<n;++i){auto score=posterior(x[i],f,gaussian);double best=*std::max_element(score.begin(),score.end());if(!std::isfinite(best)){f.labels[i]=-1;continue;}double sum=0;for(double s:score)sum+=std::exp(s-best);
        f.likelihood+=weights[i]*(best+std::log(sum));f.labels[i]=int(std::max_element(score.begin(),score.end())-score.begin());support[f.labels[i]]++;
        for(int c=0;c<k;++c)f.probability[i][c]=std::exp(score[c]-best)/sum;}
    // At least three independent observations; never a fraction of genome size.
    if(std::any_of(support.begin(),support.end(),[](size_t m){return m<3;}))return Fit{};
    // Likelihood mass and penalty both use independent-observation effective N.
    f.likelihood*=neff/n;
    double parameters=k*modeledFeatures*2+k-1;
    if(!gaussian){
        // Hard-assignment spherical Gaussian likelihood, with one shared fitted
        // variance and fitted mixing proportions. Missing cells add no evidence.
        double residual=0,cells=0,mixing=0;std::vector<double> assigned(k);
        for(size_t i=0;i<n;++i)if(f.labels[i]>=0){size_t c=size_t(f.labels[i]);assigned[c]+=weights[i];
            for(size_t j=0;j<d;++j)if(std::isfinite(x[i][j])&&std::isfinite(f.means[c][j])){
                double delta=x[i][j]-f.means[c][j];residual+=weights[i]*delta*delta;cells+=weights[i];}}
        double variance=std::max(1e-10,residual/std::max(1.,cells));
        for(size_t c=0;c<assigned.size();++c){double m=assigned[c];if(m>0)mixing+=m*std::log(m/n);f.prior[c]=m/n;
            std::fill(f.variance[c].begin(),f.variance[c].end(),variance);}
        f.likelihood=(mixing-.5*(cells*std::log(6.283185307179586*variance)+residual/variance))*neff/n;
        parameters=k*modeledFeatures+k; // means, K-1 proportions, shared variance
    }
    f.bic=-2*f.likelihood+parameters*std::log(std::max(2.,neff));
    f.aic=-2*f.likelihood+2*parameters;
    return f;
}
struct Candidate {int k=0,attempt=0;double bic=INFINITY,aic=INFINITY;bool converged=false;int iterations=0;};
inline Fit search(const Rows& x,int lo,int hi,int attempts,int iterations,uint64_t seed,bool gaussian=false,const std::vector<double>& weights={},unsigned threads=1,
                  const std::string& criterion="bic",std::vector<Candidate>* diagnostics=nullptr,unsigned maxWorkers=0,bool exhaustive=false){
    if(lo<1||hi<lo||hi>128||attempts<1||iterations<1||!threads||(criterion!="bic"&&criterion!="aic"))throw std::invalid_argument("Invalid continuous K search options");
    // Screen all K once. Keep only the winning posterior workspace between phases;
    // diagnostics contain scalars. Spend remaining restarts on the leading K.
    struct Job{int k,a;};std::vector<Job> jobs;
    for(int k=lo;k<=hi&&size_t(k)*3<=x.size();++k)jobs.push_back({k,0});
    std::vector<Candidate> records;Fit best;size_t bestJob=0;unsigned workers=std::min(threads,maxWorkers?maxWorkers:threads);
    auto phase=[&]{
        std::vector<Candidate> phaseRecords(jobs.size());std::mutex mutex;std::atomic<size_t> next{0};
        auto work=[&]{for(size_t j;(j=next.fetch_add(1))<jobs.size();){auto job=jobs[j];
            auto f=fit(x,job.k,iterations,seed+1009*job.k+104729*job.a,gaussian,weights,std::max(1u,threads/std::max(1u,std::min(workers,unsigned(jobs.size())))));
            phaseRecords[j]={job.k,job.a,f.bic,f.aic,f.converged,f.iterations};
            size_t key=size_t(job.k)*size_t(attempts)+size_t(job.a);std::lock_guard<std::mutex> lock(mutex);
            if(std::isfinite(f.score(criterion))&&(f.score(criterion)<best.score(criterion)||(f.score(criterion)==best.score(criterion)&&key<bestJob))){best=std::move(f);bestJob=key;}}};
        std::vector<std::future<void>> pool;for(size_t t=1;t<std::min<size_t>(workers,jobs.size());++t)pool.push_back(std::async(std::launch::async,work));work();for(auto& w:pool)w.get();
        records.insert(records.end(),phaseRecords.begin(),phaseRecords.end());
    };
    phase();auto ranking=records;
    std::stable_sort(ranking.begin(),ranking.end(),[&](const Candidate& a,const Candidate& b){return (criterion=="aic"?a.aic:a.bic)<(criterion=="aic"?b.aic:b.bic);});
    jobs.clear();for(size_t i=0;i<(exhaustive?ranking.size():std::min<size_t>(2,ranking.size()));++i)for(int a=1;a<attempts;++a)jobs.push_back({ranking[i].k,a});
    phase();if(diagnostics)*diagnostics=std::move(records);
    if(!std::isfinite(best.score(criterion)))throw std::runtime_error("No supported continuous fit in candidate K range");return best;
}
// Convert nonnegative sums to log rates. Equal feature scaling leaves row amplitude intact.
inline Rows features(const Rows& sums,const Rows& exposure,const std::vector<double>& sourceMass={}){
    Rows x=sums;if(x.size()!=exposure.size())throw std::invalid_argument("Exposure rows");
    for(size_t i=0;i<x.size();++i){if(x[i].size()!=exposure[i].size())throw std::invalid_argument("Exposure columns");double m=sourceMass.empty()?1:sourceMass[i];
        for(size_t j=0;j<x[i].size();++j)x[i][j]=exposure[i][j]>0&&m>0?std::log1p(sums[i][j]/(exposure[i][j]*m)):NAN;}
    for(size_t j=0;!x.empty()&&j<x[0].size();++j){double mean=0,s=0,n=0;for(const auto& row:x)if(std::isfinite(row[j])){mean+=row[j];n++;}if(n)mean/=n;for(const auto& row:x)if(std::isfinite(row[j]))s+=(row[j]-mean)*(row[j]-mean);s=n?std::sqrt(s/n):0;
        for(auto& row:x)row[j]=s>1e-10&&std::isfinite(row[j])?(row[j]-mean)/s:NAN;}return x;
}
inline void rowStandardize(Rows& x,const Rows& exposure){
    for(size_t i=0;i<x.size();++i){double mass=0,mean=0,var=0;size_t n=0;for(size_t j=0;j<x[i].size();++j)if(std::isfinite(x[i][j])&&exposure[i][j]>0){mass+=exposure[i][j];mean+=exposure[i][j]*x[i][j];++n;}if(mass>0)mean/=mass;
        for(size_t j=0;j<x[i].size();++j)if(std::isfinite(x[i][j])&&exposure[i][j]>0)var+=exposure[i][j]*(x[i][j]-mean)*(x[i][j]-mean);double sd=mass>0?std::sqrt(var/mass):0;
        for(double& v:x[i])v=n>=3&&sd>1e-6&&std::isfinite(v)?(v-mean)/sd:NAN;}
}
struct Correspondence {size_t viewA,stateA,viewB,stateB,support;double fractionA,fractionB;};
struct Reconciliation {
    std::vector<Correspondence> candidates;
    std::vector<std::vector<size_t>> mapping;
    std::vector<std::vector<size_t>> members;
    std::vector<int> parent;
    std::vector<bool> provisional;
    bool limitReached=false;
};
// Each track uses the same independent observation universe; -1 means unavailable.
// Only consistent high-overlap correspondences merge; one-sided containment is lineage.
inline Reconciliation reconcile(const std::vector<std::vector<int>>& tracks,const std::vector<double>& weights={},double threshold=.8,size_t minSupport=3,size_t limit=256){
    Reconciliation r;if(tracks.empty())return r;size_t n=tracks[0].size();std::vector<size_t> owner,local;
    for(size_t v=0;v<tracks.size();++v){if(tracks[v].size()!=n)throw std::invalid_argument("Track alignment");int k=0;for(int c:tracks[v])k=std::max(k,c+1);r.mapping.emplace_back(size_t(k));for(int c=0;c<k;++c){r.mapping.back()[c]=owner.size();owner.push_back(v);local.push_back(size_t(c));}}
    size_t m=owner.size();std::vector<double> shared(m*m),ab(m*m),ba(m*m);std::vector<size_t> common(m*m);if(!weights.empty()&&weights.size()!=n)throw std::invalid_argument("Reconciliation weights");
    for(size_t a=0;a<m;++a)for(size_t b=a+1;b<m;++b)if(owner[a]!=owner[b]){double x=0,na=0,nb=0;size_t both=0;
        for(size_t i=0;i<n;++i){int ca=tracks[owner[a]][i],cb=tracks[owner[b]][i];if(ca<0||cb<0)continue;double w=weights.empty()?1:weights[i];bool ia=size_t(ca)==local[a],ib=size_t(cb)==local[b];na+=ia?w:0;nb+=ib?w:0;if(ia&&ib){x+=w;both++;}}
        ab[a*m+b]=ba[b*m+a]=na?x/na:0;ba[a*m+b]=ab[b*m+a]=nb?x/nb:0;shared[a*m+b]=shared[b*m+a]=std::min(ab[a*m+b],ba[a*m+b]);common[a*m+b]=common[b*m+a]=both;}
    struct Edge{size_t a,b;double score;};std::vector<Edge> edges;for(size_t a=0;a<m;++a)for(size_t b=a+1;b<m;++b)if(common[a*m+b]>=minSupport&&shared[a*m+b]>=threshold)edges.push_back({a,b,shared[a*m+b]});
    std::stable_sort(edges.begin(),edges.end(),[](const Edge&a,const Edge&b){return a.score>b.score;});std::vector<size_t> group(m);std::iota(group.begin(),group.end(),0);std::vector<std::vector<size_t>> members(m);for(size_t a=0;a<m;++a)members[a]={a};
    for(const auto& e:edges){size_t a=group[e.a],b=group[e.b];if(a==b)continue;bool ok=true;
        for(size_t u:members[a])for(size_t v:members[b]){if(owner[u]==owner[v]){ok=false;continue;}bool comparable=false;for(size_t i=0;i<n&&!comparable;++i)comparable=tracks[owner[u]][i]>=0&&tracks[owner[v]][i]>=0&&(size_t(tracks[owner[u]][i])==local[u]||size_t(tracks[owner[v]][i])==local[v]);if(comparable&&(shared[u*m+v]<threshold||common[u*m+v]<minSupport))ok=false;}
        if(ok){for(size_t v:members[b]){group[v]=a;members[a].push_back(v);}members[b].clear();}}
    std::vector<size_t> remap(m);for(size_t a=0;a<m;++a)if(!members[a].empty()){size_t g=r.members.size();for(size_t v:members[a])remap[v]=g;r.members.push_back(members[a]);}
    r.parent.assign(r.members.size(),-1);r.provisional.assign(r.members.size(),true);r.limitReached=r.members.size()>limit;
    for(size_t g=0;g<r.members.size();++g){r.provisional[g]=r.members[g].size()<2;double best=0;for(size_t u:r.members[g])for(size_t v=0;v<m;++v)if(remap[v]!=g&&common[u*m+v]>=minSupport&&ab[u*m+v]>=threshold&&ba[u*m+v]<threshold&&ba[u*m+v]>best){best=ba[u*m+v];r.parent[g]=int(remap[v]);}}
    for(size_t g=0;g<r.parent.size();++g){std::vector<bool> seen(r.parent.size());int p=int(g);while(p>=0&&!seen[size_t(p)]){seen[size_t(p)]=true;p=r.parent[size_t(p)];}if(p>=0)r.parent[g]=-1;}
    for(size_t a=0;a<m;++a)for(size_t b=a+1;b<m;++b)if(common[a*m+b])r.candidates.push_back({owner[a],local[a],owner[b],local[b],common[a*m+b],ab[a*m+b],ba[a*m+b]});
    for(auto& map:r.mapping)for(size_t& g:map)g=remap[g];return r;
}
inline void refineScores(std::vector<double>& score,const Reconciliation& r){
    std::vector<size_t> order(score.size());std::iota(order.begin(),order.end(),0);
    auto depth=[&](size_t g){size_t d=0;for(int p=r.parent[g];p>=0&&d<r.parent.size();p=r.parent[size_t(p)])++d;return d;};
    std::stable_sort(order.begin(),order.end(),[&](size_t a,size_t b){return depth(a)<depth(b);});
    for(size_t p:order)if(score[p]>0){double children=0;for(size_t g=0;g<score.size();++g)if(r.parent[g]==int(p)&&score[g]>0)children+=score[g];
        if(children>1e-6&&children>=.25*score[p]){double inherited=score[p];for(size_t g=0;g<score.size();++g)if(r.parent[g]==int(p)&&score[g]>0)score[g]+=inherited*score[g]/children;score[p]=0;}}
}
inline std::vector<int> consensus(const std::vector<std::vector<int>>& tracks,const Reconciliation& r,const std::vector<double>& viewWeight={}){
    if(tracks.empty())return {};std::vector<int> result(tracks[0].size(),-1);
    for(size_t i=0;i<result.size();++i){std::vector<double> score(r.members.size());for(size_t v=0;v<tracks.size();++v)if(tracks[v][i]>=0)score[r.mapping[v][tracks[v][i]]]+=viewWeight.empty()?1:viewWeight[v];
        refineScores(score,r);
        if(!score.empty()){double best=*std::max_element(score.begin(),score.end());size_t ties=0;for(double s:score)ties+=best>0&&std::abs(s-best)<1e-12;if(ties==1)result[i]=int(std::max_element(score.begin(),score.end())-score.begin());
            else if(ties>1){for(size_t v=0;v<tracks.size();++v)if(tracks[v][i]>=0){size_t g=r.mapping[v][tracks[v][i]];if(std::abs(score[g]-best)<1e-12){result[i]=int(g);break;}}}}}return result;
}
inline bool unresolved(int a,int b,const Reconciliation& r){
    if(a<0||b<0)return false;size_t x=r.mapping[0][size_t(a)],y=r.mapping[1][size_t(b)];return x!=y&&r.parent[x]!=int(y)&&r.parent[y]!=int(x);
}
inline void lineage(const std::string& path,const Reconciliation& r){std::ofstream f(path);if(!f)throw std::runtime_error("Cannot write lineage");f<<"state\tparent\tprovisional\tsupporting_fits\n";for(size_t g=0;g<r.members.size();++g)f<<g+1<<'\t'<<r.parent[g]+1<<'\t'<<r.provisional[g]<<'\t'<<r.members[g].size()<<'\n';
    auto slash=path.find_last_of("/\\");std::ofstream c((slash==std::string::npos?std::string{}:path.substr(0,slash+1))+"state_correspondence.tsv");
    c<<"view_a\tstate_a\tview_b\tstate_b\tindependent_shared_observations\tfraction_a\tfraction_b\trelation\n";
    for(const auto& edge:r.candidates){size_t a=r.mapping[edge.viewA][edge.stateA],b=r.mapping[edge.viewB][edge.stateB];c<<edge.viewA+1<<'\t'<<edge.stateA+1<<'\t'<<edge.viewB+1<<'\t'<<edge.stateB+1<<'\t'<<edge.support<<'\t'<<edge.fractionA<<'\t'<<edge.fractionB<<'\t'<<(a==b?"matched":r.parent[a]==int(b)||r.parent[b]==int(a)?"refinement":"unresolved")<<'\n';}
}
}
