//
// Created by Yuanyuan Qi on 8/23/23.
//

#ifndef BBT_BBT_H
#define BBT_BBT_H

#include <list>
#include <vector>

class Input;
class AncestryGraph;

typedef std::list<std::pair<int,int> > BBT;
//void bulid_AD_pairs(int n, BBT &A);

bool BBT_check(const Input &input, const BBT &bbt);
bool BBT_check_AD(const Input &input, const BBT &bbt );

struct BBT_hc{
    std::vector<std::list<BBT> > hc;
    std::vector<int> order_of_mutations;
    std::vector<std::vector<bool> > binary_m;
    const AncestryGraph & GF;
    explicit BBT_hc(const AncestryGraph &GF);
    int _max_depth,_max_count;
    void _hc(int depth);
    void h_c(int limit=0x7fffffff);
};
#endif //BBT_BBT_H
