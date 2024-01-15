//
// Created by Yuanyuan Qi on 7/6/23.
//

#include "Input.h"

#include <fstream>
#include <list>

Input::Input(const std::string &filename, int r): r(r) {
    std::ifstream fin(filename);
    fin >> m >> n;
    data = std::vector<std::vector<double> >(m, std::vector<double>(n) );
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; j++) {
            fin >> data[i][j];
        }
    }
    update_r();
}

Input::Input(int m, int n, int r): m(m),n(n),r(r), data(m, std::vector<double>(n) ) {
}

void Input::update_r() {
    if (r < 0){
        std::vector<double> sof(n,0);
        for (int i = 0; i < m; i++){
            for (int j = 0; j < n; j++){
                sof[j]+=data[i][j];
            }
        }
        double max_f = -1;
        for (int i=0; i<n; i++){
            if (sof[i] > max_f ){
                max_f = sof[i];
                this->r = i;
            }
        }
    }
}

Input::Input(const std::vector<std::vector<double> > & f,int r):m(f.size()),n(f[0].size()),r(r){
    data=f;
}

//void Input::Rsof() {
//    sof.resize(n, 0);
//    for (int i = 0; i < m; i++) {
//        for (int j = 0; j < n; j++) {
//            sof[j] += data[i][j];
//        }
//    }
//    if (r >= 0) return;
//    double maxf = -1e300;
//    for (int j = 0; j < n; j++) {
//        if (maxf < sof[j]) {
//            maxf = sof[j];
//            r = j;
//        }
//    }
//}


Reads::Reads(const std::string &var_file, const std::string &ref_file, const std::string &info_file) {
    std::ifstream f1(var_file);
    std::ifstream f2(ref_file);
    std::ifstream f0(info_file);
    f0 >> m >> n_mutations >> n;
    cluster = std::vector<int>(n_mutations);
    cnt_cluster = std::vector<int>(n, 0);
    var = std::vector<std::vector<int> >(m, std::vector<int>(n_mutations));
    ref = var;
    tot = var;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n_mutations; ++j) {
            f1 >> var[i][j];
            f2 >> ref[i][j];
            tot[i][j] = var[i][j] + ref[i][j];
        }
    }
    for (int i = 0; i < n_mutations; i++) {
        f0 >> cluster[i];
        cnt_cluster[cluster[i]]++;
    }
    f0 >> r;

    std::vector<std::vector<int> > va(m, std::vector<int>(n, 0)), re(m, std::vector<int>(n, 0));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n_mutations; j++) {
            va[i][cluster[j]] += var[i][j];
            re[i][cluster[j]] += ref[i][j];
        }
    }
    summary_freq = std::vector<std::vector<double> > (m, std::vector<double>(n, 0));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            summary_freq[i][j] = (double) va[i][j] / (va[i][j] + re[i][j]);
        }
    }

    std::vector<std::vector<std::vector<int> > > cl_tot_list
            (m, std::vector<std::vector<int> >(n, std::vector<int>(n_mutations + 1, 0)));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n_mutations; j++) {
            int &index(cl_tot_list[i][cluster[j]][0]);
            cl_tot_list[i][cluster[j]][++index] = tot[i][j];
        }
    }

    tot_c_med = std::vector<std::vector<int> > (m, std::vector<int>(n, 0));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            std::sort(cl_tot_list[i][j].begin() + 1, cl_tot_list[i][j].begin() + 1 + cl_tot_list[i][j][0]);
            int med_index_l = cl_tot_list[i][j][0] / 2 + 1;
            int med_index_r = (cl_tot_list[i][j][0] + 1) / 2;
            tot_c_med[i][j] = cl_tot_list[i][j][med_index_l] + cl_tot_list[i][j][med_index_r];
            tot_c_med[i][j] /= 2;
        }
    }
}

Input* Reads::Bootstrap(std::mt19937 & rng) const {
    auto *p = new Input(m, n, r);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            std::binomial_distribution distr(tot_c_med[i][j],summary_freq[i][j]);
            int k = distr(rng);
            p->data[i][j] = (double)k/tot_c_med[i][j];
        }
    }
    p->update_r();
    return p;
}

//Input *Reads::All() const {
//    auto *p = new Input(m, n, r);
//    std::vector<std::vector<int> > va(m, std::vector<int>(n, 0)), re(m, std::vector<int>(n, 0));
//    for (int i = 0; i < m; i++) {
//        for (int j = 0; j < n_mutations; j++) {
//            va[i][cluster[j]] += var[i][j];
//            re[i][cluster[j]] += ref[i][j];
//        }
//    }
//    for (int i = 0; i < m; i++) {
//        for (int j = 0; j < n; j++) {
//            p->data[i][j]=(double)va[i][j]/(va[i][j]+re[i][j]);
//        }
//    }
//    p->update_r();
//    return p;
//}
