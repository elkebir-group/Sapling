//
// Created by Yuanyuan Qi on 6/14/23.
//

#ifndef UNIPPM_ANCESTRYGRAPH_H
#define UNIPPM_ANCESTRYGRAPH_H

#include <vector>
#include <list>

class Input;

struct AncestryGraph {
    const Input & data;
    explicit AncestryGraph(const Input &data, bool init=true);
    std::vector<std::list<int> > possible_children;
    std::vector<std::list<int> > possible_parent;
    std::vector<std::vector<int> > arc_set_index;
    std::vector<std::pair<int,int> > arc_set;
};

typedef std::list<std::pair<int,int> > BBT;

struct AncestryGraph_BBT: public AncestryGraph{
    std::vector<std::vector<double > >  remainder;
    std::vector<int> parent;
    AncestryGraph_BBT(const Input &data,  const BBT &bbt);
};

#endif //UNIPPM_ANCESTRYGRAPH_H
