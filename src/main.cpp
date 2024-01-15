//
// Created by Yuanyuan Qi on 8/25/23.
//

#include "Input.h"
#include "AncestryGraph.h"
#include "BBT.h"

float EPS;

int main(int argc, char ** argv) {
    EPS = 1e-4;
    Reads *praw;
    int seed = 1;
    if (argc >= 5) seed = std::stoi(argv[4]);
    std::mt19937 rng(seed);
    praw = new Reads(argv[2], argv[3], argv[1]);
    auto *p = praw->Bootstrap(rng);
    Input data(*p);
//    Input data("input.txt",0);
    AncestryGraph GF(data);
    BBT_hc main_hc(GF);
    main_hc.h_c();
    for (auto &a: main_hc.hc) {
        printf("%d\n", a.size());
        for (auto &T:a){
            printf("[");
            for(auto arc:T.arcs)
                printf("(%d,%d),",arc.first,arc.second);
            printf("]\n");
        }
    }
    delete praw;
    delete p;
    return 0;
}