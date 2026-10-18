#pragma once
// Shared native trans preprocessing. Keep the copy in Mixer include/slice in sync.
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>
namespace hic_trans {
struct Balance {
    std::vector<double> rowFactors, columnFactors;
    int iterations=0;
    double relativeError=0, meanContact=0;
};
inline void memoryGuard(size_t rows,size_t columns,double gb,double copies=12) {
    if(!std::isfinite(gb)||gb<=0)throw std::invalid_argument("Memory budget must be positive");
    if(double(rows)*double(columns)*copies*8>gb*1073741824.)
        throw std::runtime_error("MEMORY_BUDGET: matrix workspace exceeds --max-memory-gb");
}
// exposure[i][j] is eligible target mass, before summing source rows into pools.
// It is additive when pools are formed. rowMass selects equal-pool or genomic-mass targets.
template<class Data> Balance balance(Data& data,int maxIterations=1000,double tolerance=1e-4) {
    size_t n=data.values.size(),p=data.weights.size();
    if(!n||!p||maxIterations<1||!std::isfinite(tolerance)||tolerance<=0)throw std::invalid_argument("Invalid balancing dimensions or tolerance");
    if(data.exposure.empty()) {
        data.exposure.assign(n,typename decltype(data.exposure)::value_type(p,0));
        for(size_t i=0;i<n;++i)for(size_t j=0;j<p;++j)if(std::isfinite(data.values[i][j]))data.exposure[i][j]=data.weights[j];
    }
    if(data.rowMass.empty())data.rowMass.assign(n,1);
    if(data.exposure.size()!=n||data.rowMass.size()!=n||(!data.sourceSize.empty()&&data.sourceSize.size()!=n))throw std::invalid_argument("Invalid exposure dimensions");
    std::vector<double> rt(n),ct(p),rs(n),cs(p);
    for(size_t i=0;i<n;++i){
        if(data.values[i].size()!=p||data.exposure[i].size()!=p||!(data.rowMass[i]>0)||!std::isfinite(data.rowMass[i])||(!data.sourceSize.empty()&&(!(data.sourceSize[i]>0)||!std::isfinite(data.sourceSize[i]))))throw std::invalid_argument("Invalid balancing row");
        for(size_t j=0;j<p;++j)if(std::isfinite(data.values[i][j])) {
            double x=data.values[i][j],e=data.exposure[i][j];
            if(x<0||!std::isfinite(e)||e<0)throw std::invalid_argument("Invalid trans count or exposure");
            if(e>0&&data.weights[j]>0){rs[i]+=x;cs[j]+=x;}
        }
    }
    Balance out;out.rowFactors.assign(n,1);out.columnFactors.assign(p,1);
    double total=0,opportunity=0;
    for(size_t i=0;i<n;++i)for(size_t j=0;j<p;++j){
        auto& x=data.values[i][j];double e=data.exposure[i][j];
        if(!std::isfinite(x)||!(e>0)||!(rs[i]>0)||!(cs[j]>0)||data.weights[j]<=0){x=std::numeric_limits<float>::quiet_NaN();data.exposure[i][j]=0;continue;}
        // Remove additive pool size from per-source exposure, then apply chosen mass.
        // sourceSize is always the represented count, independent of weighting policy.
        double sourceSize=data.sourceSize.empty()?1:data.sourceSize[i];
        double target=e/sourceSize*data.rowMass[i];
        rt[i]+=target;ct[j]+=target;total+=x;opportunity+=target;
    }
    if(!(total>0)||!(opportunity>0))throw std::runtime_error("No usable trans contacts for balancing");
    out.meanContact=total/opportunity;
    for(size_t i=0;i<n;++i)if(!(rt[i]>0))out.rowFactors[i]=0;
    for(size_t j=0;j<p;++j)if(!(ct[j]>0))out.columnFactors[j]=0;
    for(int it=1;it<=maxIterations;++it){
        std::fill(cs.begin(),cs.end(),0);
        for(size_t i=0;i<n;++i)if(rt[i]>0){double sum=0;for(size_t j=0;j<p;++j)if(std::isfinite(data.values[i][j]))sum+=data.values[i][j]*out.columnFactors[j];
            if(!(sum>0)||!std::isfinite(sum))throw std::runtime_error("Invalid row marginal during trans balancing");
            out.rowFactors[i]=rt[i]/sum;
            for(size_t j=0;j<p;++j)if(std::isfinite(data.values[i][j]))cs[j]+=data.values[i][j]*out.rowFactors[i];}
        for(size_t j=0;j<p;++j)if(ct[j]>0){if(!(cs[j]>0)||!std::isfinite(cs[j]))throw std::runtime_error("Invalid column marginal during trans balancing");out.columnFactors[j]=ct[j]/cs[j];}
        double error=0;
        for(size_t i=0;i<n;++i)if(rt[i]>0){double sum=0;for(size_t j=0;j<p;++j)if(std::isfinite(data.values[i][j]))sum+=data.values[i][j]*out.columnFactors[j];error=std::max(error,std::abs(sum*out.rowFactors[i]/rt[i]-1));}
        out.iterations=it;out.relativeError=error;
        if(error<=tolerance)break;
        if(it==maxIterations)throw std::runtime_error("Compressed trans balancing did not converge; residual "+std::to_string(error));
    }
    for(size_t i=0;i<n;++i)for(size_t j=0;j<p;++j)if(std::isfinite(data.values[i][j])){
        // Express balanced features per original source bin. Cleaning divides by
        // per-source eligible target mass, so unequal pools have comparable rates.
        double value=data.values[i][j]*out.rowFactors[i]*out.columnFactors[j]*out.meanContact/data.rowMass[i];
        if(!std::isfinite(value)||value>std::numeric_limits<float>::max())throw std::runtime_error("Balanced count exceeds float range");
        data.values[i][j]=float(value);
        double sourceSize=data.sourceSize.empty()?1:data.sourceSize[i];
        data.exposure[i][j]/=sourceSize;
    }
    return out;
}
// Log1p rates retain measured zeros. Winsorization limits influence without
// turning biological extremes into missing observations. Column z, optional
// weighted row z, then sqrt(feature weight), in that order.
template<class Data> void clean(Data& data,bool rowZ=false,bool clip=true) {
    auto& x=data.values;size_t n=x.size(),p=data.weights.size();
    if(!n||!p)throw std::runtime_error("Empty continuous feature matrix");
    if(data.intervals.size()!=n||(!data.exposure.empty()&&data.exposure.size()!=n)||(!data.rowMass.empty()&&data.rowMass.size()!=n))throw std::invalid_argument("Invalid continuous-feature metadata");
    for(size_t i=0;i<n;++i)if(x[i].size()!=p||(!data.exposure.empty()&&data.exposure[i].size()!=p))throw std::invalid_argument("Ragged continuous-feature matrix");
    std::vector<double> masses=data.rowMass;if(masses.empty())masses.assign(n,1);
    for(size_t i=0;i<n;++i)for(size_t j=0;j<p;++j){double e=data.exposure.empty()?data.weights[j]:data.exposure[i][j];
        if(!std::isfinite(x[i][j])||x[i][j]<0||!(e>0))x[i][j]=std::numeric_limits<float>::quiet_NaN();
        else x[i][j]=float(std::log1p(double(x[i][j])/e));}
    for(size_t j=0;j<p;++j){std::vector<float> values;for(const auto& row:x)if(std::isfinite(row[j]))values.push_back(row[j]);
        if(values.empty())continue;
        float lo=0,hi=0;
        if(clip){size_t h=size_t(.99*(values.size()-1));lo=0;std::nth_element(values.begin(),values.begin()+h,values.end());hi=values[h];size_t middle=values.size()/2;std::nth_element(values.begin(),values.begin()+middle,values.end());hi=std::max(hi,float(values[middle]+std::log(100.)));}
        double sum=0,mass=0,var=0;
        for(size_t i=0;i<n;++i)if(std::isfinite(x[i][j])){if(clip)x[i][j]=std::clamp(x[i][j],lo,hi);sum+=masses[i]*x[i][j];mass+=masses[i];}
        double mean=sum/std::max(1.,mass);
        for(size_t i=0;i<n;++i)if(std::isfinite(x[i][j]))var+=masses[i]*(x[i][j]-mean)*(x[i][j]-mean);
        double sd=std::sqrt(var/std::max(1.,mass));
        for(auto& row:x)if(std::isfinite(row[j]))row[j]=sd>1e-8?float((row[j]-mean)/sd):std::numeric_limits<float>::quiet_NaN();
    }
    if(rowZ)for(size_t i=0;i<n;++i){auto& row=x[i];double sum=0,mass=0,var=0;size_t features=0;
        for(size_t j=0;j<p;++j)if(std::isfinite(row[j])&&data.weights[j]>0){double w=data.exposure.empty()?data.weights[j]:data.exposure[i][j];sum+=w*row[j];mass+=w;++features;}
        double mean=sum/std::max(1.,mass);
        for(size_t j=0;j<p;++j)if(std::isfinite(row[j]))var+=(data.exposure.empty()?data.weights[j]:data.exposure[i][j])*(row[j]-mean)*(row[j]-mean);
        double sd=std::sqrt(var/std::max(1.,mass));
        for(auto& value:row)if(std::isfinite(value))value=features>=3&&sd>1e-6?float((value-mean)/sd):std::numeric_limits<float>::quiet_NaN();}
    size_t kept=0;
    for(size_t i=0;i<n;++i){size_t good=0;for(float v:x[i])good+=std::isfinite(v);if(!good)continue;
        for(size_t j=0;j<p;++j)x[i][j]*=std::sqrt(data.weights[j]);
        if(kept!=i){x[kept]=std::move(x[i]);data.intervals[kept]=data.intervals[i];if(!data.exposure.empty())data.exposure[kept]=std::move(data.exposure[i]);if(!data.rowMass.empty())data.rowMass[kept]=data.rowMass[i];if(!data.sourceSize.empty())data.sourceSize[kept]=data.sourceSize[i];}++kept;}
    x.resize(kept);data.intervals.resize(kept);if(!data.exposure.empty())data.exposure.resize(kept);if(!data.rowMass.empty())data.rowMass.resize(kept);if(!data.sourceSize.empty())data.sourceSize.resize(kept);
    if(!kept)throw std::runtime_error("No rows have variable features after normalization");
}
}
