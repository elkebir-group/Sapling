//
// Created by Yuanyuan Qi on 8/23/23.
//

#include <queue>

#include "BBT.h"
#include "Input.h"
#include "AncestryGraph.h"

#include "gurobi_c++.h"

extern float EPS;

bool BBT_check(const Input& data, const BBT &bbt){
    GRBEnv env(true);
    env.set(GRB_IntParam_OutputFlag, 0);
    env.start();
    GRBModel model(env);
    AncestryGraph_BBT GF(data,bbt);
    std::vector<GRBVar> arc(GF.arc_set.size());
    model.set(GRB_IntParam_MIPFocus, GRB_MIPFOCUS_FEASIBILITY);
    model.set(GRB_IntParam_LazyConstraints, 1);
    for (int i = 0; i < GF.arc_set.size(); i++) {
        arc[i] = model.addVar(0, 1, 0, GRB_INTEGER,
                              "arc[" + std::to_string(GF.arc_set[i].first) + "," +
                              std::to_string(GF.arc_set[i].second)+"]");
    }

    GRBLinExpr sum;
    for (int i = 0; i < data.m; i++) {
        for (int p = 0; p < data.n; p++) {
            for (auto q: GF.possible_children[p]) {
                sum += arc[GF.arc_set_index[p][q]]*data.data[i][q];
            }
            model.addConstr(GF.remainder[i][p]+EPS >= sum );
            sum.clear();
        }
    }

    //////one incoming edge
    for (int i = 0; i < data.n; i++) {
        if(i==data.r or GF.parent[i]>=0) continue;
        for (auto j: GF.possible_parent[i]) {
            sum += arc[GF.arc_set_index[j][i]];
        }
        model.addConstr(sum == 1);
        sum.clear();
    }

    model.optimize();
    if (model.get(GRB_IntAttr_Status)==GRB_INFEASIBLE)
        return false;
    return true;
}

//void bulid_AD_pairs(int n, BBT &A){
//    if (!A.ad_pairs.empty()) return;
//    std::vector<std::list<int> > link_list(n);
//    std::queue<std::pair<int,int> > Q;
//    for (auto it:A.arcs){
//        link_list[it.first].push_back(it.second);
//        Q.push(it);
//    }
//    std::pair<int,int> fr;
//    while (!Q.empty()){
//        fr = Q.front();
//        A.ad_pairs.push_back(fr);
//        Q.pop();
//        for (auto ne :link_list[fr.second]){
//            Q.emplace(std::pair(fr.first,ne));
//        }
//    }
//}

bool BBT_check_AD(const Input& data, const BBT &bbt){
    GRBEnv env(true);
    env.set(GRB_IntParam_OutputFlag, 0);
    env.start();
    GRBModel model(env);
    AncestryGraph GF(data);
    std::vector<GRBVar> arc(GF.arc_set.size());
    model.set(GRB_IntParam_MIPFocus, GRB_MIPFOCUS_FEASIBILITY);
    model.set(GRB_IntParam_LazyConstraints, 1);
    for (int i = 0; i < GF.arc_set.size(); i++) {
        arc[i] = model.addVar(0, 1, 0, GRB_INTEGER,
                              "arc[" + std::to_string(GF.arc_set[i].first) + "," +
                              std::to_string(GF.arc_set[i].second)+"]");
    }

    GRBLinExpr sum;
    for (int i = 0; i < data.m; i++) {
        for (int p = 0; p < data.n; p++) {
            for (auto q: GF.possible_children[p]) {
                sum += arc[GF.arc_set_index[p][q]]*data.data[i][q];
            }
            model.addConstr(data.data[i][p]+EPS >= sum );
            sum.clear();
        }
    }

    //////one incoming edge
    for (int i = 0; i < data.n; i++) {
        if(i==data.r) continue;
        for (auto j: GF.possible_parent[i]) {
            sum += arc[GF.arc_set_index[j][i]];
        }
        model.addConstr(sum == 1);
        sum.clear();
    }

    ////////AD constraint
    std::vector<std::vector<GRBVar> > ad(data.n,std::vector<GRBVar>(data.n));
    std::vector<std::vector<std::vector<GRBVar> > >
            adp(data.n,std::vector<std::vector<GRBVar> >(data.n,std::vector<GRBVar>(data.n)));

    for (int i = 0; i < data.n; i++){
        for (int j = 0; j < data.n; j++){
            if (i==j) continue;
            ad[i][j] = model.addVar(0,1,0,GRB_CONTINUOUS,
                                          "ad["+std::to_string(i)+","+std::to_string(j)+"]");
            for (int mid = 0; mid < data.n; mid++){
                if (i==mid || mid==j) continue;
                adp[i][j][mid] = model.addVar(0,1,0,GRB_CONTINUOUS,
                                              "ad["+std::to_string(i)+","+std::to_string(j)
                                              +"] passing_last("+std::to_string(mid)+")");
            }
        }
    }

    for (int i = 0; i < data.n; i++){
        for (int j = 0; j < data.n; j++){
            if (i==j) continue;
            if (GF.arc_set_index[i][j]>=0) {
                sum += arc[GF.arc_set_index[i][j]];
                model.addConstr(ad[i][j] >= arc[GF.arc_set_index[i][j]]);
            }
            for (int mid = 0; mid < data.n; mid++){
                if (i==mid || mid==j) continue;
                if (GF.arc_set_index[mid][j] < 0 ) continue;
                sum+=adp[i][j][mid];
                model.addConstr(ad[i][j] >= adp[i][j][mid]);
                model.addConstr(adp[i][j][mid]>=ad[i][mid]+arc[GF.arc_set_index[mid][j]]-1);
                model.addConstr(adp[i][j][mid]<=ad[i][mid]);
                model.addConstr(adp[i][j][mid]<=arc[GF.arc_set_index[mid][j]]-1);
            }
            model.addConstr(ad[i][j]<=sum);
            sum.clear();
        }
    }

//    bulid_AD_pairs(data.n,bbt);

    for (auto it:bbt.arcs){
        model.addConstr(ad[it.first][it.second]>=1);
    }

    model.optimize();
    if (model.get(GRB_IntAttr_Status)==GRB_INFEASIBLE)
        return false;
    return true;
}

BBT_hc::BBT_hc(const AncestryGraph &GF): hc(GF.data.n), GF(GF), order_of_mutations(GF.data.n),
binary_m(GF.data.n, std::vector<bool>(GF.data.n,false)) {
    order_of_mutations[0] = this->GF.data.r;
    hc[0].push_back(BBT());
    binary_m[0][GF.data.r] = true;
}

void BBT_hc::_hc(int depth) {
    std::list<int> front;
    bool flag;
    for (int i = 0; i < GF.data.n; i++) {
        if (binary_m[depth][i]) continue;
        flag = true;
        for (auto j: GF.possible_parent[i]) {
            if (!binary_m[depth][j]) {
                flag = false;
                break;
            }
        }
        if (flag) front.push_back(i);
    }
    int opt = -1;
    std::list<BBT> res;
    BBT tmp;
    for (auto i: front) {
        for (auto &t: hc[depth]) {
            for (auto p: GF.possible_parent[i]) {
                tmp = t;
                tmp.arcs.emplace_back(p, i);
                if (BBT_check(GF.data, tmp))
                    res.push_back(tmp);
            }
        }
        if (opt<0 or hc[depth+1].size()>res.size()){
            opt = i;
            hc[depth+1]=res;
            res.clear();
        }
    }
    order_of_mutations[depth+1]=opt;
    std::copy(binary_m[depth].begin(),binary_m[depth].end(),binary_m[depth+1].begin());
    binary_m[depth+1][opt]=true;
}

void BBT_hc::h_c(int limit) {
    for (int i=0; i<GF.data.n-1; i++){
        _hc(i);
        if (hc[i+1].size()>=limit) break;
    }
}
