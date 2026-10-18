#include "slice/slice.hpp"
#include <cmath>
#include <stdexcept>
#include <iostream>
void check(bool b) { if (!b) throw std::runtime_error("Numerical regression"); }
int main() {
    slice::Random r(0); check(r.nextLong()==uint64_t(-4962768465676381896LL));
    slice::Row a{1,2,slice::nan,4},b{2,4,3,8};
    check(std::abs(slice::squaredDistance(a,b)-28.)<1e-6);
    check(std::abs(slice::absoluteDistance(a,b)-28./3)<1e-6);
    check(std::abs(slice::correlation(a,b)-1)<1e-6);
    slice::Matrix data{{0,1},{0,2},{0,3},{0,100}};
    auto means=slice::kmeans(data,1,200,42,1,false);
    auto medians=slice::kmeans(data,1,200,42,2,true);
    check(means[0].center[1]==26.5f); check(medians[0].center[1]==2.5f);
    check(std::isnan(slice::correlation({1,1},{2,2})));
    std::cout<<"numerical tests passed\n";
}
