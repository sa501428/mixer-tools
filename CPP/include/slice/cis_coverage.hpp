#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <ostream>
#include <string>
#include <vector>

namespace cis_coverage {
constexpr double cutoff = -2.5;
struct Profile {
    std::vector<double> nonzero, z;
    std::vector<bool> bad;
    double mean=0, sd=0;
};
// Population moments over nonempty rows, independently per chromosome.
// Empty rows fail even when the positive-row distribution has zero variance.
inline Profile evaluate(std::vector<double> nonzero) {
    Profile p;p.nonzero=std::move(nonzero);p.z.resize(p.nonzero.size());p.bad.resize(p.nonzero.size());
    size_t n=0;double m2=0;
    for(double value:p.nonzero)if(value>0){++n;double delta=value-p.mean;p.mean+=delta/n;m2+=delta*(value-p.mean);}
    p.sd=n?std::sqrt(std::max(0.,m2/n)):0;
    for(size_t i=0;i<p.nonzero.size();++i){
        p.z[i]=p.sd>0?(p.nonzero[i]-p.mean)/p.sd:0;
        p.bad[i]=p.nonzero[i]<=0||p.z[i]<cutoff;
    }
    return p;
}
inline void reportHeader(std::ostream& out) {
    out<<"chromosome\tstart\tend\tcis_nonzero\tzscore\tmean_nonzero\tsd_nonzero\tretained\n";
}
inline void report(std::ostream& out,const std::string& chrom,int64_t length,int resolution,const Profile& p) {
    for(size_t i=0;i<p.nonzero.size()&&int64_t(i)*resolution<length;++i)
        out<<chrom<<'\t'<<int64_t(i)*resolution<<'\t'<<std::min(length,int64_t(i+1)*resolution)<<'\t'
           <<p.nonzero[i]<<'\t'<<p.z[i]<<'\t'<<p.mean<<'\t'<<p.sd<<'\t'<<!p.bad[i]<<'\n';
}
} // namespace cis_coverage
