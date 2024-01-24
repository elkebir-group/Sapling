//
// Created by Yuanyuan Qi on 1/11/24.
//

//#include "Input.h"
//#include "AncestryGraph.h"
//#include "BBT.h"
#include "LLH.h"
#include <vector>
#include <list>
#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>

float EPS;

namespace py = pybind11;

//struct Fcase{
//    Input * F;
////    AncestryGraph * G;
//    explicit Fcase(const py::list & F);
//    bool verify_BBT(const py::list & edges);
//};
//
//Fcase::Fcase(const py::list & F) {
//    std::vector<std::vector<double> > casted(F.size());
//    int idx=0;
//    for (auto & i:F){
//        casted[idx++] = i.cast<std::vector<double> >();
//    }
//    this->F = new Input(casted);
//    this->F->update_r();
////    G = new AncestryGraph(*this->F);
//}
//
//bool Fcase::verify_BBT(const py::list & edges){
//    Tree casted;
//    for(auto & edge:edges){
//        casted.push_back(edge.cast<std::pair<int,int> >());
//    }
//    auto a = BBT_check_AD(*F,casted);
//    return a;
//}

struct maxLLH{
    LLH *llh;
    std::vector<std::vector<double> > *rp_usage, *rp_f_snv;
    Tree *T;
    maxLLH(const py::list& Var, const py::list& Tot, int n_intervals, double alpha);
    double LLH_SNV_(const py::list & T);
    double LLH_SNV_t_(const py::list & partial_T);
    ~maxLLH();
};

maxLLH::maxLLH(const py::list & Var, const py::list & Tot, int n_intervals, double alpha)
{
    IM casted_var(Var.size()),casted_tot(Tot.size());
    int idx=0;
    for(auto & var_r:Var){
        casted_var[idx++] = var_r.cast<std::vector<int> >();
    }
    idx=0;
    for(auto & tot_r:Tot){
        casted_tot[idx++] = tot_r.cast<std::vector<int> >();
    }
    llh = new LLH(casted_var,casted_tot,n_intervals,alpha);
    rp_usage = new std::vector<std::vector<double> >();
    rp_f_snv = new std::vector<std::vector<double> >();
    T = new Tree ();
}

double maxLLH::LLH_SNV_(const py::list & T){

    Tree casted;
    std::vector<int> tmp;
    for(auto & edge:T){
        casted.push_back(edge.cast<std::pair<int,int> >());
    }
    return llh->LLH_SNV(casted,*rp_f_snv,*rp_usage);
}

double maxLLH::LLH_SNV_t_(const py::list &partial_T) {
    Tree casted;
    std::vector<int> tmp;
    for(auto & edge:partial_T){
        casted.push_back(edge.cast<std::pair<int,int> >());
    }
    return llh->LLH_SNV_t(casted,*rp_f_snv, *rp_usage, *T);
}

maxLLH::~maxLLH() {
    delete llh;
    delete rp_usage;
    delete rp_f_snv;
    delete T;
}


PYBIND11_MODULE(pyBBT, MODULE) {
//    py::class_<Fcase>(MODULE, "case").def(py::init<py::list >())
//            .def("verify_BBT", &Fcase::verify_BBT);
    py::class_<maxLLH>(MODULE, "maxllh").def(py::init<py::list,py::list,int,double>())
            .def("llh_snv",&maxLLH::LLH_SNV_);
}