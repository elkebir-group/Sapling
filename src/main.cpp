//
// Created by Yuanyuan Qi on 8/25/23.
//

#include "Input.h"
#include "AncestryGraph.h"
#include "BBT.h"
#include "LLH.h"

float EPS;
//
// Created by Yuanyuan Qi on 1/11/24.
//
int main(int argc, char ** argv) {
    EPS = 1e-4;
    Reads *praw;
    int seed = 1;
    if (argc >= 5) seed = std::stoi(argv[4]);
    std::mt19937 rng(seed);
    praw = new Reads(argv[2], argv[3], argv[1]);
    LLH llh(praw->var,praw->tot,100,1);
//    Input data(*p);
////    Input data("input.txt",0);
//    AncestryGraph GF(data);
//    BBT_hc main_hc(GF);
//    main_hc.h_c();
//    for (auto &a: main_hc.hc) {
//        printf("%d\n", a.size());
//        for (auto &T:a){
//            printf("[");
//            for(auto arc:T)
//                printf("(%d,%d),",arc.first,arc.second);
//            printf("]\n");
//        }
//    }
//    delete p;
    std::vector<std::vector<double> > result_usage,result_f;
    Tree T = {{0,2},{0,3}},res_T;
//    llh.LLH_SNV(T,result_usage,result_f);
    llh.LLH_SNV_adt(T,result_usage,result_f,res_T);
    for (auto a:result_f){
        for(auto f:a){
            printf("%.3lf ",f);
        }
        printf("\n");
    }
    printf("Tree:[");
    for (auto a:res_T){
        printf("(%d,%d)",a.first,a.second);
    }
    printf("]\n");
    delete praw;
    return 0;
}