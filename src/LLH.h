//
// Created by Yuanyuan Qi on 1/23/21.
//

#ifndef BBT_LLH_H
#define BBT_LLH_H

#include <vector>
#include <list>

typedef std::list<std::pair<int,int> > Tree;
typedef std::vector<std::vector<int> > IM;

class LLH {
public:
    std::vector<std::vector<double> > F_upper, F_lower;
    std::vector<std::vector<std::vector<double> > > split;
    IM var,ref,tot;
    Tree opt;
    int m,n;
    std::vector<std::list<int> > PP,PC;

    int n_intervals;
    int arc_size;
    std::vector<std::pair<int,int> > arc_set;
    std::vector<std::vector<int> > arc_set_idx;

//    double alpha_SNV,alpha_Meth;

//    static void generate_splits(int var, int ref, double lower, double uppper, double alpha, int n_interval, std::vector<double> & split);

    LLH(IM Var, IM Tot, int n_intervals=100, double alpha=1);//, double alpha_snv=-1, double alpha_meth=-1);
    double LLH_SNV(const Tree & T,
                   std::vector<std::vector<double> > & result_usage,
                   std::vector<std::vector<double> > & result_f_snv);

    double LLH_SNV_t(const Tree & T,
                   std::vector<std::vector<double> > & result_usage,
                   std::vector<std::vector<double> > & result_f_snv,
                   Tree & result_t);
    double LLH_SNV_adt(const Tree & T,
                     std::vector<std::vector<double> > & result_usage,
                     std::vector<std::vector<double> > & result_f_snv,
                     Tree & result_t);
};


#endif //BBT_LLH_H
