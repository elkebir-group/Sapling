//
// Created by Yuanyuan Qi on 6/14/23.
//

#include "AncestryGraph.h"
#include "Input.h"

extern float EPS;

AncestryGraph::AncestryGraph(const Input &data, bool init):
data(data),possible_children(data.n), possible_parent(data.n),
arc_set_index(data.n,std::vector<int>(data.n,-1)) {
    if (!init) return;
    int count = 0;
    bool tmp_flag;
    for (auto i = 0; i < data.n; ++i) {
        for (auto j = 0; j < data.n; ++j) {
            if (i == j or j == data.r) continue;
            tmp_flag = true;
            for (auto _sample = 0; _sample < data.m; _sample++) {
                if (data.data[_sample][i] + EPS < data.data[_sample][j]) { //i can't be ancestor of j
                    tmp_flag = false;
                    break;
                }
            }
            if (tmp_flag) {
                count += 1;
                possible_children[i].push_back(j);
                possible_parent[j].push_back(i);
            }
        }
    }
    arc_set.resize(count);
    for (auto i = 0, _count = 0; i < data.n; ++i) {
        for (auto j: possible_children[i]) {
            arc_set_index[i][j] = _count;
            arc_set[_count] = {i, j};
            ++_count;
        }
    }
    for (int i = 0; i < data.n; i++) {
        if (i != data.r && possible_parent[i].empty()){
            throw std::string("Invalid root in Ancestry Graph.");
        }
    }
}

AncestryGraph_BBT::AncestryGraph_BBT(const Input &data, const BBT &bbt):
AncestryGraph(data, false),
remainder(data.m,std::vector<double>(data.n)),
parent(data.n,-1) {
    int count = 0;
    bool tmp_flag;
    for (auto i = 0; i < data.m; ++i) {
        for (auto j = 0; j < data.n; ++j) {
            remainder[i][j] = data.data[i][j];
        }
    }
    for (auto a:bbt){
        parent[a.second]=a.first;
        for (auto i = 0; i < data.m; ++i) {
            remainder[i][a.first]-=data.data[i][a.second];
        }
    }
    for (auto i = 0; i < data.n; ++i) {
        for (auto j = 0; j < data.n; ++j) {
            if (i == j or j==data.r or parent[j]>=0) continue;
            tmp_flag = true;
            for (auto _sample = 0; _sample < data.m; _sample++) {
                if (remainder[_sample][i]+EPS < data.data[_sample][j]) { //i can't be ancestor of j
                    tmp_flag = false;
                    break;
                }
            }
            if (tmp_flag) {
                count += 1;
                possible_children[i].push_back(j);
                possible_parent[j].push_back(i);
            }
        }
    }
    arc_set.resize(count);
    for (auto i = 0, _count = 0; i < data.n; ++i) {
        for (auto j: possible_children[i]) {
            arc_set_index[i][j] = _count;
            arc_set[_count] = {i, j};
            ++_count;
        }
    }
}
