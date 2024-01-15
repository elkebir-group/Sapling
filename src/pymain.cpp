//
// Created by Yuanyuan Qi on 1/11/24.
//

#include "Input.h"
#include "AncestryGraph.h"
#include "BBT.h"
#include <pybind11/pybind11.h>

float EPS;

namespace py = pybind11;

struct Fcase{
    Input * F;
//    AncestryGraph * G;
    Fcase(int seed, const std::vector<std::vector<double> > & F);
    bool verify_BBT(const std::list<std::pair<int,int> > & edges);
};

Fcase::Fcase(int seed, const std::vector<std::vector<double> > &F) {
    this->F = new Input(F);
    this->F->update_r();
//    G = new AncestryGraph(*this->F);
}

bool Fcase::verify_BBT(const std::list<std::pair<int, int> > &edges){
    auto a = BBT_check_AD(*F,BBT({edges}));
}

PYBIND11_MODULE(pyBBT, m) {
    py::class_<Fcase>(m,"case").def(py::init<int,const std::vector<std::vector<double> > >())
                        .def("verify_BBT", &Fcase::verify_BBT);
}