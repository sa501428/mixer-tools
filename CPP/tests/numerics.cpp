#include "slice/slice.hpp"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <random>
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
    // Compare bounded distances with the exhaustive implementation, including
    // chromosome-shaped missingness, ties and centroid support changes.
    std::mt19937_64 rng(42);std::normal_distribution<float> noise(0,1);
    for(int seed=0;seed<12;++seed){slice::Matrix matrix(240,slice::Row(36));
        for(size_t i=0;i<matrix.size();++i)for(size_t j=0;j<matrix[i].size();++j)
            matrix[i][j]=(j/6==i%6)?slice::nan:float(i%4)*3+noise(rng);
        matrix.back().assign(36,slice::nan);
        auto exhaustive=slice::kmeans(matrix,4,80,seed,1,true,false);
        auto bounded=slice::kmeans(matrix,4,80,seed,4,true,true);
        check(exhaustive.size()==bounded.size());
        for(size_t c=0;c<exhaustive.size();++c){check(exhaustive[c].members==bounded[c].members);
            for(size_t j=0;j<36;++j)check(exhaustive[c].center[j]==bounded[c].center[j]||(std::isnan(exhaustive[c].center[j])&&std::isnan(bounded[c].center[j])));}}
    std::cout<<"numerical tests and exact median bounds passed\n";
}
