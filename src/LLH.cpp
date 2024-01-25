//
// Created by Yuanyuan Qi on 1/23/21.
//

#include "LLH.h"
#include <gurobi_c++.h>

#include <boost/math/special_functions/beta.hpp>
#include <stack>

#define EPS 1e-4

inline double log_(double x,double eps){
    if(x<eps){
        return (log(eps)-(eps-x)*eps*eps);
    }
    return log(x);
}

inline std::pair<double,double> interval_from_reads(int var,int ref, double alpha, bool fix_01=false, double n_fix=0.5) {
    double lower=0, upper=1;
    if (alpha >= 1) {
        return {0, 1};
    }
    lower = boost::math::ibeta_inv(var + n_fix, ref + n_fix, (1 - alpha) / 2);
    upper = boost::math::ibeta_inv(var + n_fix, ref + n_fix, (1 + alpha) / 2);
    if (fix_01) {
        if (var == 0) {
            lower = 0;
        }
        if (ref == 0) {
            upper = 1;
        }
    }
    return {lower, upper};
}

LLH::LLH(IM Var, IM Tot, int n_intervals, double alpha):
m(Var.size()),
n(Var[0].size()),
n_intervals(n_intervals),
var((Var)),
tot((Tot)),
ref(Var.size(),std::vector<int>(Var[0].size())),
F_upper(Var.size(),std::vector<double>(Var[0].size(),0)),
F_lower(Var.size(),std::vector<double>(Var[0].size(),0)),
PP(Var[0].size()),
PC(Var[0].size()),
arc_set_idx(Var[0].size(),std::vector<int>(Var[0].size(),-1)),
arc_size(0)
{
    split.resize(m,
                 std::vector<std::vector<double> >(n, std::vector<double>(n_intervals + 1)));

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            ref[i][j] = tot[i][j]-var[i][j];
            auto tmp = interval_from_reads(var[i][j],ref[i][j],alpha);
            F_upper[i][j] = tmp.second;
            F_lower[i][j] = tmp.first;
        }
    }

    for (int p = 0; p < n; p++) {
        for (int q = 0; q < n; q++) {
            if (p == q) continue;
            bool flag = true;
            for (int i = 0; i < m; i++) {
                if (F_lower[i][q]>F_upper[i][p]){
                    flag = false;
                    break;
                }
            }
            if (flag){
                PP[q].push_back(p);
                PC[p].push_back(q);
                arc_set_idx[p][q]=arc_size++;
            }
        }
    }

    arc_set.resize(arc_size);
    for (int p = 0; p < n; p++) {
        for(auto q:PC[p]){
            arc_set[arc_set_idx[p][q]]={p,q};
        }
    }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n_intervals + 1; k++) {
                split[i][j][k] = F_lower[i][j] + (F_upper[i][j] - F_lower[i][j]) / (n_intervals) * k;
            }
        }
    }
}

double LLH::LLH_SNV(const Tree &T, std::vector<std::vector<double> > & result_usage,
                    std::vector<std::vector<double> > & result_f_snv) {
    result_f_snv.resize(m,std::vector<double>(n));
    result_usage.resize(m,std::vector<double>(n));
    GRBEnv env(true);
//    env.set(GRB_IntParam_OutputFlag, 0);
    env.start();
    GRBModel model(env);
    std::vector<std::vector<GRBVar> >
            f(m,std::vector<GRBVar>(n));
    std::vector<std::vector<std::vector<GRBVar> > >
            lambda(m,std::vector<std::vector<GRBVar> >
            (n,std::vector<GRBVar> (n_intervals+1) ));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            f[i][j] = model.addVar(F_lower[i][j], F_upper[i][j], 0, GRB_CONTINUOUS,
                                   "f["+std::to_string(i)+"]["+std::to_string(j)+"]");
            for(int k = 0; k <= n_intervals; k++){
                lambda[i][j][k] = model.addVar(0, 1, 0, GRB_CONTINUOUS,
                                               "lambda["+std::to_string(i)+"]["+
                                               std::to_string(j)+"]["+
                                               std::to_string(k)+"]");
            }
        }
    }

    std::vector<std::list<int> > TPC(n);
    for (auto e: T){
        TPC[e.first].push_back(e.second);
    }

    GRBLinExpr sum;
    for (int i = 0; i < m; i++){
        for (int j = 0 ; j < n; j++){
            for (auto item:TPC[j])
                sum += f[i][item];
            model.addConstr(f[i][j] >= sum);
            sum.clear();
        }
    }

    GRBLinExpr sum2;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k <= n_intervals; k++){
                sum += lambda[i][j][k];
                sum2 += lambda[i][j][k] * split[i][j][k];
            }
            model.addConstr(sum == 1);
            model.addConstr(sum2 == f[i][j]);
            sum.clear();
            sum2.clear();
        }
    }

    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            for (int k = 0; k <= n_intervals; k++){
                sum += lambda[i][j][k] * (var[i][j]*log_(split[i][j][k],EPS)+
                                          ref[i][j]*log_(1-split[i][j][k],EPS));
            }
        }
    }

    model.setObjective(sum, GRB_MAXIMIZE);
    model.optimize();

    auto obj_val_inferred = model.get(GRB_DoubleAttr_ObjVal);
//    double log_binom_coeffs = 0;
//    double obj_val = 0;
//    if (debug){
//        for (int i = 0; i < m; i++) {
//            for (int j = 0; j < n; j++) {
//                log_binom_coeffs += logcomb(var[i][j] + ref[i][j],
//                                            var[i][j]);
//                obj_val += logbinom(var[i][j] + ref[i][j],
//                                    var[i][j],
//                                    f[i][j].get(GRB_DoubleAttr_X), EPS);
//            }
//        }
////        std::cout << "GUROBI: " << obj_val_inferred + log_binom_coeffs << " -- recomputed: " << obj_val << std::endl;
//    }


    ////get freq and result_usage
    for(int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            result_f_snv[i][j] = f[i][j].get(GRB_DoubleAttr_X);
        }
    }

    for(int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            result_usage[i][j] = result_f_snv[i][j];
            for (auto ch: TPC[j]){
                result_usage[i][j]-=result_f_snv[i][ch];
            }
        }
    }

    return obj_val_inferred;
}




class opt_T_callback: public GRBCallback{
public:
    LLH * master_class;
    std::vector<GRBVar> & arc_vars;
    opt_T_callback(LLH & llh, std::vector<GRBVar> & arc_v):master_class(&llh),arc_vars(arc_v){}
    void callback() {
        if (where != GRB_CB_MIPSOL) return;
        std::list<std::pair<int,int> > edges,cycle;
        std::vector<bool> instack(master_class->n,false),visited(master_class->n,false);
        std::stack<std::pair<int,std::list<int>::iterator> > st;
        int idx=0;
        for (auto e:master_class->arc_set){
            auto val = getSolution(arc_vars[idx]);
            idx++;
            if (val>0.99){
                edges.push_back(e);
            }
        }

        std::vector<std::list<int> > TPC(master_class->n);
        for (auto e:edges){
            TPC[e.first].push_back(e.second);
        }
        bool flag=false;
        int next,loop_p;
        for (int i=0;i<master_class->n;i++){
            if (!visited[i]){
                st.push({i,TPC[i].begin()});
                instack[i]=true;
                visited[i]=true;
            }
            while (!st.empty()){
                if(st.top().second==TPC[st.top().first].end()){
                    instack[st.top().first] = false;
                    st.pop();
                    continue;
                }
                next = *st.top().second;
                st.top().second++;
                if (instack[next]){
                    flag=true;
                    break;
                }
                st.push({next,TPC[next].begin()});
                instack[next] = true;
                visited[next] = true;
            }
            if (flag) break;
        }
        if (flag) {
            loop_p = next;
            while (st.top().first != loop_p) {
                cycle.emplace_back(st.top().first, next);
                next = st.top().first;
                st.pop();
            }
            cycle.emplace_back(st.top().first, next);
            GRBLinExpr sum;
            for (auto e: cycle) {
                sum += arc_vars[master_class->arc_set_idx[e.first][e.second]];
            }
            addLazy(sum <= cycle.size() - 1);
        }
    }
};

double LLH::LLH_SNV_t(const Tree & T,
                 std::vector<std::vector<double> > & result_usage,
                 std::vector<std::vector<double> > & result_f_snv,
                 Tree & result_t) {
    result_f_snv.resize(m,std::vector<double>(n));
    result_usage.resize(m,std::vector<double>(n));
    result_t.clear();
    GRBEnv env(true);
//    env.set(GRB_IntParam_OutputFlag, 0);
    env.start();
    GRBModel model(env);
    std::vector<std::vector<GRBVar> > f(m, std::vector<GRBVar>(n));
    std::vector<std::vector<GRBVar> > ef(m, std::vector<GRBVar>(arc_size));
    std::vector<std::vector<std::vector<GRBVar> > >
            lambda(m, std::vector<std::vector<GRBVar> >
            (n, std::vector<GRBVar>(n_intervals + 1)));
    std::vector<GRBVar> isroot(n);
    std::vector<GRBVar> arc(arc_size);

    model.set(GRB_IntParam_LazyConstraints, 1);

    //////f vars and splits
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            f[i][j] = model.addVar(F_lower[i][j], F_upper[i][j], 0, GRB_CONTINUOUS,
                                   "f["+std::to_string(i)+"]["+std::to_string(j)+"]");
            for(int k = 0; k <= n_intervals; k++){
                lambda[i][j][k] = model.addVar(0, 1, 0, GRB_CONTINUOUS,
                                               "lambda["+std::to_string(i)+"]["+
                                               std::to_string(j)+"]["+
                                               std::to_string(k)+"]");
            }
        }
    }

    GRBLinExpr sum,sum2;

    //////pw split of f
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k <= n_intervals; k++){
                sum += lambda[i][j][k];
                sum2 += lambda[i][j][k] * split[i][j][k];
            }
            model.addConstr(sum == 1);
            model.addConstr(sum2 == f[i][j]);
            sum.clear();
            sum2.clear();
        }
    }

    //////edge vars
    for (int i = 0; i < arc_size; i++) {
        arc[i] = model.addVar(0, 1, 0, GRB_INTEGER,
                              "arc[" + std::to_string(arc_set[i].first) + "," +
                              std::to_string(arc_set[i].second) + "]");
    }

    //////edge vars * f vars
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < arc_size; j++) {
            ef[i][j] = model.addVar(0,1,0, GRB_CONTINUOUS,
                                    "f["+std::to_string(i)+","+
                                    std::to_string(arc_set[j].second)+"]*arc["+
                                    std::to_string(arc_set[j].first) + "," +
                                    std::to_string(arc_set[j].second) + "]");
            model.addConstr(ef[i][j] >= arc[j]+f[i][arc_set[j].second]-1);
            model.addConstr(ef[i][j] <= f[i][arc_set[j].second]);
            model.addConstr(ef[i][j] <= arc[j]);
        }
    }

    ///////sum condition
    for (int i = 0; i < m; i++) {
        for (int p = 0; p < n; p++) {
            for (auto q: PC[p]) {
                sum += ef[i][arc_set_idx[p][q]];
            }
            model.addConstr(f[i][p] >= sum);
            sum.clear();
        }
    }

    //////one root
    for (int p = 0; p < n; p++) {
        isroot[p] = model.addVar(0,1,0,GRB_INTEGER,
                                 "r["+std::to_string(p)+"]");
        sum+=isroot[p];
    }
    model.addConstr(sum==1);
    sum.clear();

    //////one incoming edge
    for (int i = 0; i < n; i++) {
        for (auto j: PP[i]) {
            sum += arc[arc_set_idx[j][i]];
        }
        model.addConstr(sum+isroot[i] == 1);
        sum.clear();
    }


    //////objective value
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            for (int k = 0; k <= n_intervals; k++){
                sum += lambda[i][j][k] * (var[i][j]*log_(split[i][j][k],EPS)+
                                          ref[i][j]*log_(1-split[i][j][k],EPS));
            }
        }
    }

    //////cycle prevention
    for (auto e:T){
        model.addConstr(arc[arc_set_idx[e.first][e.second]]==1);
    }
    model.setObjective(sum, GRB_MAXIMIZE);
    opt_T_callback callback(*this,arc);
    model.setCallback(&callback);
    model.optimize();
//    model.computeIIS();
//    model.write("iis.ilp");

    auto obj_val_inferred = model.get(GRB_DoubleAttr_ObjVal);

    for(int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            result_f_snv[i][j] = f[i][j].get(GRB_DoubleAttr_X);
        }
    }

    for(int i=0;i<arc_size;i++){
        if (arc[i].get(GRB_DoubleAttr_X)>0.99){
            result_t.push_back(arc_set[i]);
        }
    }


    std::vector<std::list<int> > TPC(n);
    for (auto e: result_t){
        TPC[e.first].push_back(e.second);
    }

    for(int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            result_usage[i][j] = result_f_snv[i][j];
            for (auto ch: TPC[j]){
                result_usage[i][j]-=result_f_snv[i][ch];
            }
        }
    }

    return obj_val_inferred;
}



double LLH::LLH_SNV_adt(const Tree & T,
                      std::vector<std::vector<double> > & result_usage,
                      std::vector<std::vector<double> > & result_f_snv,
                      Tree & result_t) {
    result_f_snv.resize(m,std::vector<double>(n));
    result_usage.resize(m,std::vector<double>(n));
    result_t.clear();
    GRBEnv env(true);
//    env.set(GRB_IntParam_OutputFlag, 0);
    env.start();
    GRBModel model(env);
    std::vector<std::vector<GRBVar> > f(m, std::vector<GRBVar>(n));
    std::vector<std::vector<GRBVar> > ef(m, std::vector<GRBVar>(arc_size));
    std::vector<std::vector<std::vector<GRBVar> > >
            lambda(m, std::vector<std::vector<GRBVar> >
            (n, std::vector<GRBVar>(n_intervals + 1)));
    std::vector<GRBVar> isroot(n);
    std::vector<GRBVar> arc(arc_size);

//    model.set(GRB_IntParam_LazyConstraints, 1);

    //////f vars and splits
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            f[i][j] = model.addVar(F_lower[i][j], F_upper[i][j], 0, GRB_CONTINUOUS,
                                   "f["+std::to_string(i)+"]["+std::to_string(j)+"]");
            for(int k = 0; k <= n_intervals; k++){
                lambda[i][j][k] = model.addVar(0, 1, 0, GRB_CONTINUOUS,
                                               "lambda["+std::to_string(i)+"]["+
                                               std::to_string(j)+"]["+
                                               std::to_string(k)+"]");
            }
        }
    }

    GRBLinExpr sum,sum2;

    //////pw split of f
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k <= n_intervals; k++){
                sum += lambda[i][j][k];
                sum2 += lambda[i][j][k] * split[i][j][k];
            }
            model.addConstr(sum == 1);
            model.addConstr(sum2 == f[i][j]);
            sum.clear();
            sum2.clear();
        }
    }

    //////edge vars
    for (int i = 0; i < arc_size; i++) {
        arc[i] = model.addVar(0, 1, 0, GRB_INTEGER,
                              "arc[" + std::to_string(arc_set[i].first) + "," +
                              std::to_string(arc_set[i].second) + "]");
    }

    //////edge vars * f vars
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < arc_size; j++) {
            ef[i][j] = model.addVar(0,1,0, GRB_CONTINUOUS,
                                    "f["+std::to_string(i)+","+
                                    std::to_string(arc_set[j].second)+"]*arc["+
                                    std::to_string(arc_set[j].first) + "," +
                                    std::to_string(arc_set[j].second) + "]");
            model.addConstr(ef[i][j] >= arc[j]+f[i][arc_set[j].second]-1);
            model.addConstr(ef[i][j] <= f[i][arc_set[j].second]);
            model.addConstr(ef[i][j] <= arc[j]);
        }
    }

    ///////sum condition
    for (int i = 0; i < m; i++) {
        for (int p = 0; p < n; p++) {
            for (auto q: PC[p]) {
                sum += ef[i][arc_set_idx[p][q]];
            }
            model.addConstr(f[i][p] >= sum);
            sum.clear();
        }
    }

    //////one root
    for (int p = 0; p < n; p++) {
        isroot[p] = model.addVar(0,1,0,GRB_INTEGER,
                                 "r["+std::to_string(p)+"]");
        sum+=isroot[p];
    }
    model.addConstr(sum==1);
    sum.clear();

    //////one incoming edge
    for (int i = 0; i < n; i++) {
        for (auto j: PP[i]) {
            sum += arc[arc_set_idx[j][i]];
        }
        model.addConstr(sum+isroot[i] == 1);
        sum.clear();
    }

    //////TODO: needs changing
    std::vector<std::vector<GRBVar> > ad(n,std::vector<GRBVar>(n));
    std::vector<std::vector<std::vector<GRBVar> > >
            adp(n,std::vector<std::vector<GRBVar> >(n,std::vector<GRBVar>(n)));

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i==j) continue;
            ad[i][j] = model.addVar(0,1,0,GRB_CONTINUOUS,
                                    "ad["+std::to_string(i)+","+std::to_string(j)+"]");
            for (int mid = 0; mid < n; mid++){
                if (i==mid || mid==j) continue;
                adp[i][j][mid] = model.addVar(0,1,0,GRB_CONTINUOUS,
                                              "ad["+std::to_string(i)+","+std::to_string(j)
                                              +"]_passing_last("+std::to_string(mid)+")");
            }
        }
    }

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i==j) continue;
            if (arc_set_idx[i][j]>=0) {
                sum += arc[arc_set_idx[i][j]];
                model.addConstr(ad[i][j] >= arc[arc_set_idx[i][j]]);
            }
            for (int mid = 0; mid < n; mid++){
                if (i==mid || mid==j) continue;
                if (arc_set_idx[mid][j] < 0 ) continue;
                sum+=adp[i][j][mid];
                model.addConstr(ad[i][j] >= adp[i][j][mid]);
                model.addConstr(adp[i][j][mid]>=ad[i][mid]+arc[arc_set_idx[mid][j]]-1);
                model.addConstr(adp[i][j][mid]<=ad[i][mid]);
                model.addConstr(adp[i][j][mid]<=arc[arc_set_idx[mid][j]]);
            }
            model.addConstr(ad[i][j]<=sum);
            sum.clear();
        }
    }


    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
            model.addConstr(ad[i][j]+ad[j][i]<=1);
        }
    }

    for (auto e:T){
        model.addConstr(ad[e.first][e.second]==1);
    }


    //////objective value
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            for (int k = 0; k <= n_intervals; k++){
                sum += lambda[i][j][k] * (var[i][j]*log_(split[i][j][k],EPS)+
                                          ref[i][j]*log_(1-split[i][j][k],EPS));
            }
        }
    }

    model.setObjective(sum, GRB_MAXIMIZE);
//    opt_T_callback callback(*this,arc);
//    model.setCallback(&callback);
    model.optimize();

    auto obj_val_inferred = model.get(GRB_DoubleAttr_ObjVal);

    for(int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            result_f_snv[i][j] = f[i][j].get(GRB_DoubleAttr_X);
        }
    }

    for(int i=0;i<arc_size;i++){
        if (arc[i].get(GRB_DoubleAttr_X)>0.99){
            result_t.push_back(arc_set[i]);
        }
    }

    std::vector<std::list<int> > TPC(n);
    for (auto e: result_t){
        TPC[e.first].push_back(e.second);
    }

    for(int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            result_usage[i][j] = result_f_snv[i][j];
            for (auto ch: TPC[j]){
                result_usage[i][j]-=result_f_snv[i][ch];
            }
        }
    }

    return obj_val_inferred;
}
