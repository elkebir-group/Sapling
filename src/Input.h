//
// Created by Yuanyuan Qi on 7/6/23.
//

#ifndef UNIPPM_INPUT_H
#define UNIPPM_INPUT_H

#include <vector>
#include <string>
#include <random>

struct Input {
    int n,m; //n mutations, m samples
    int r;
    std::vector<std::vector<double> >  data;
    explicit Input(const std::string& filename, int r = -1);
    explicit Input(const std::vector<std::vector<double> > & f, int r=-1);
    Input(int n, int m, int r = -1);
    void update_r();
//    std::vector<double> sof;
//    void Rsof();
};

struct Reads {
    int n_mutations;
    int n,m;//n_clusters
    int r;
    std::vector<std::vector<int> > var,ref,tot;
    std::vector<std::vector<int> > tot_c_med;
    std::vector<std::vector<double> > summary_freq;
    std::vector<int> cluster,cnt_cluster;
    Reads(const std::string& var_file,const std::string & ref_file, const std::string & info_file);
    Input* Bootstrap(std::mt19937 & rng) const;
//    Input* All() const;
};
#endif //UNIPPM_INPUT_H
