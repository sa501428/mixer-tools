#include "slice/cis_coverage.hpp"
#include <iostream>
#include <stdexcept>
int main(){
    auto check=[](bool ok){if(!ok)throw std::runtime_error("cis coverage contract failed");};
    std::vector<double> rows(100,100);rows[0]=1;rows[1]=0;
    auto p=cis_coverage::evaluate(rows);
    check(p.bad[0]&&p.z[0]<-2.5&&p.bad[1]);
    for(size_t i=2;i<rows.size();++i)check(!p.bad[i]);
    check(std::abs(p.mean-99)<1e-10&&std::abs(p.sd-std::sqrt(98.))<1e-10);
    // Empty gaps must not distort the positive-row population.
    rows.insert(rows.end(),1000,0);auto gaps=cis_coverage::evaluate(rows);
    check(gaps.mean==p.mean&&gaps.sd==p.sd);
    auto flat=cis_coverage::evaluate({10,10,0,10});
    check(flat.sd==0&&!flat.bad[0]&&!flat.bad[1]&&flat.bad[2]&&!flat.bad[3]);
    auto empty=cis_coverage::evaluate({0,0});check(empty.bad[0]&&empty.bad[1]);
    check(cis_coverage::evaluate({}).bad.empty());
    std::cout<<"Cis nonzero z-score, empty rows, gaps and zero variance passed\n";
}
