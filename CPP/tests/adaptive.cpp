#include "adaptive.hpp"
#include <iostream>
#include <random>
#include <set>
using namespace hic_adaptive;
void require(bool b,const char* message){if(!b)throw std::runtime_error(message);}
int main(){
    std::vector<std::vector<int>> tracks(2,std::vector<int>(180));
    for(int i=0;i<180;++i){tracks[0][i]=i/20;tracks[1][i]=i/20;}
    // Eight broad classes, with one class split into two independently supported children.
    for(int i=0;i<180;++i)tracks[0][i]=std::min(i/20,7);
    auto r=reconcile(tracks);auto calls=consensus(tracks,r);std::set<int> states(calls.begin(),calls.end());
    require(!states.count(-1)&&states.size()==9,"A supported ninth subtype survives broad/fine reconciliation");
    require(std::count_if(r.parent.begin(),r.parent.end(),[](int p){return p>=0;})==2,"Refinement lineage has two children");
    for(int i=0;i<180;++i)if(tracks[1][i]==8)tracks[0][i]=-1;
    r=reconcile(tracks);calls=consensus(tracks,r);require(calls.back()>=0&&r.provisional[size_t(calls.back())],"View-specific supported state remains provisional");
    // Disjoint orientations can be linked by a third partition, without forcing contact-enrichment matches.
    std::vector<std::vector<int>> partial(3,std::vector<int>(40,-1));
    for(int i=0;i<40;++i){partial[i<20?0:1][i]=i%2;partial[2][i]=i%2;}
    r=reconcile(partial);require(r.members.size()==2,"Overlapping partitions resolve disjoint sides");
    // A high-overlap chain must not override contradictory direct evidence.
    std::vector<std::vector<int>> bridge(3,std::vector<int>(100,0));
    for(int i=0;i<100;++i){bridge[0][i]=i<50?0:1;bridge[1][i]=i<42?0:(i<92?1:0);bridge[2][i]=i<34?0:(i<84?1:0);}
    r=reconcile(bridge,{},.8,3);require(r.members.size()>2,"Incompatible endpoints block a transitive merge");
    auto conflict=reconcile(bridge,{},.8,3);auto retained=consensus(bridge,conflict);
    require(std::none_of(retained.begin(),retained.end(),[](int c){return c<0;}),"Unresolved identities preserve primary calls and coverage");
    std::mt19937_64 rng(47);std::normal_distribution<double> noise(0,.15);Rows data;
    for(int i=0;i<180;++i){std::vector<double> row(6);int c=i<80?0:(i<160?1:2);for(int j=0;j<6;++j)row[j]=c*4+noise(rng);if(i%5==0)row[1]=NAN;data.push_back(row);}
    auto fit=search(data,2,4,2,70,42,true);require(fit.means.size()==3,"Regularized masked GMM finds a rare third state");
    size_t correct=0;std::set<int> ids;for(int c=0;c<3;++c){std::vector<int> n(3);for(int i=0;i<180;++i)if((i<80?0:(i<160?1:2))==c)n[size_t(fit.labels[i])]++;ids.insert(int(std::max_element(n.begin(),n.end())-n.begin()));correct+=*std::max_element(n.begin(),n.end());}
    require(ids.size()==3&&correct>175,"Masked GMM recovers planted memberships");
    Rows larger;for(int r=0;r<10;++r)larger.insert(larger.end(),data.begin(),data.end());
    auto serial=hic_adaptive::fit(larger,3,50,47,true,{},1),parallel=hic_adaptive::fit(larger,3,50,47,true,{},4);
    require(serial.bic==parallel.bic&&serial.means==parallel.means&&serial.probability==parallel.probability,"Parallel GMM inference preserves deterministic row-order reductions");
    Rows same(30,std::vector<double>(4,1));require(!std::isfinite(hic_adaptive::fit(same,3,20,42,true).bic),"Collapsed components are rejected after bounded repair");
    std::cout<<"Adaptive refinement, partial matching, contradictory chains, rare states and missing-coordinate GMM passed\n";
}
