//
// Created by Yuanyuan Qi on 1/23/21.
//

#include "LLH.h"
#include <gurobi_c++.h>

#include <boost/math/special_functions/beta.hpp>
#include <utility>

#define EPS 1e-6

inline double log_(double x,double eps){
    if(x<eps){
        return (log(eps)-(eps-x)*eps*eps);
    }
    return log(x);
}

inline std::pair<double,double> interval_from_reads(int var,int ref, double alpha, bool fix_01=false, double n_fix=0.5) {
    double lower, upper;
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

LLH::LLH(IM Var, IM  Tot, int n_intervals, double alpha):
m(Var.size()),
n(Var[0].size()),
var(std::move(Var)),
tot(std::move(Tot)),
ref(Var.size(),std::vector<int>(Var[0].size())),
F_upper(Var.size(),std::vector<double>(Var[0].size(),0)),
F_lower(Var.size(),std::vector<double>(Var[0].size(),0)),
PP(Var[0].size()),
PC(Var[0].size())
{
    this->n_intervals = n_intervals;

    split.resize(m,
                 std::vector<std::vector<double> >(n, std::vector<double>(n_intervals + 1)));

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            ref[i][j] = tot[i][j]-var[i][j];
            auto tmp = interval_from_reads(var[i][j],ref[i][j],alpha);
            F_upper[i][j] = tmp.first;
            F_lower[i][j] = tmp.second;
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
            }
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
    GRBEnv env(true);
    env.set(GRB_IntParam_OutputFlag, 0);
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

    GRBLinExpr sum;
    for (int i = 0; i < m; i++){
        for (int j = 0 ; j < n; j++){
            for (auto item:PC[j])
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
    double log_binom_coeffs = 0;
    double obj_val = 0;
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
            for (auto ch: PC[j]){
                result_usage[i][j]-=result_f_snv[i][ch];
            }
        }
    }

    return obj_val;
}
