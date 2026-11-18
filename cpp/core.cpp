// Sapling-CN scoring core: the loss of many clone trees at once, on a pool of threads.
//
// A tree is given by its parent array, the SNV clusters on the edge into every clone, and
// the copy-number events on the edge into every clone ((segment, x, y): segment s acquires
// state (x, y)). Scoring a tree is the same as sapling.cn.grow.score:
//   1. the state of every clone on every segment (the nearest event at or above it);
//   2. for every cluster c on the edge into w and every segment s with SNVs of c, the
//      multiplicities (over the clones) that an SNV of c on s may have (snv_options);
//   3. alternately, every SNV takes its options' best fit to its reads in all samples given
//      the usages, and fastppm-CN fits the usages of every sample given the multiplicities
//      (started from the last usages), at most `rounds` times or until no multiplicity changes.
// The loss of a tree is the sum over the samples of the fastppm-CN loss and of the squared
// error of the calls of states that no clone has.
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <limits>
#include <mutex>
#include <set>
#include <thread>
#include <tuple>
#include <vector>

#include "Solvers/CopyNumberADMM/Solver.h"

namespace py = pybind11;
using State = std::pair<int, int>;
using Call = std::tuple<int, int, double>;  // (x, y, fraction)

namespace {

constexpr double INF = std::numeric_limits<double>::infinity();

int allele(const State& st, int h) { return h == 0 ? st.first : st.second; }

struct Tree {
    std::vector<int> parent;
    std::vector<std::vector<int>> clusters;                   // clusters on the edge into v
    std::vector<std::vector<std::tuple<int, int, int>>> cn;   // (s, x, y) on the edge into v
};

// Mutated copies in a child with xc copies of the allele, given c mutated of xp in its parent.
std::vector<int> counts(int xp, int xc, int c) {
    if (c == 0 || xp == xc) return {c};
    if (c == xp) return {xc};
    std::vector<int> out;
    if (xc > xp) {
        for (int x = c; x <= c + xc - xp; ++x) out.push_back(x);
    } else {
        for (int x = std::max(0, xc - (xp - c)); x <= std::min(c, xc); ++x) out.push_back(x);
    }
    return out;
}

using Assign = std::vector<std::pair<int, int>>;  // (clone, mutated copies)

// Port of sapling.cn.tree.snv_options (same options, same order, same cap).
struct Options {
    const std::vector<State>& st;                 // state of every clone on the segment
    const std::vector<std::vector<int>>& ch;
    int cap;

    // assignments for the strict subtree of x, given c mutated copies of allele h at x
    std::vector<Assign> below(int h, int x, int c) const {
        std::vector<std::vector<Assign>> per_child;
        for (int y : ch[x]) {
            std::vector<Assign> opts;
            for (int cy : counts(allele(st[x], h), allele(st[y], h), c)) {
                for (const Assign& rest : below(h, y, cy)) {
                    Assign a{{y, cy}};
                    a.insert(a.end(), rest.begin(), rest.end());
                    opts.push_back(std::move(a));
                    if (static_cast<int>(opts.size()) >= cap) break;
                }
            }
            per_child.push_back(std::move(opts));
        }
        std::vector<Assign> res;
        for (const auto& o : per_child)
            if (o.empty()) return {Assign{}};
        std::vector<size_t> idx(per_child.size(), 0);
        while (true) {  // cartesian product, the last factor fastest (itertools.product)
            Assign a;
            for (size_t j = 0; j < per_child.size(); ++j) a.insert(a.end(), per_child[j][idx[j]].begin(), per_child[j][idx[j]].end());
            res.push_back(std::move(a));
            if (static_cast<int>(res.size()) >= cap) break;
            int j = static_cast<int>(per_child.size()) - 1;
            while (j >= 0 && ++idx[j] == per_child[j].size()) idx[j--] = 0;
            if (j < 0) break;
        }
        return res;
    }

    // multiplicity vectors (over the clones) of an SNV introduced on the edge into w
    std::vector<std::vector<int>> of(const std::vector<int>& parent, int w, bool event_at_w) const {
        const int k = parent.size(), p = parent[w];
        std::set<std::vector<int>> vecs;
        for (int h = 0; h < 2; ++h) {
            std::vector<int> starts;
            if (allele(st[w], h) >= 1) starts.push_back(1);                       // after the change
            if (event_at_w && allele(st[p], h) >= 1)                               // before it
                for (int c : counts(allele(st[p], h), allele(st[w], h), 1))
                    if (c >= 1) starts.push_back(c);
            std::sort(starts.begin(), starts.end());
            starts.erase(std::unique(starts.begin(), starts.end()), starts.end());
            for (int cw : starts) {
                for (const Assign& rest : below(h, w, cw)) {
                    std::vector<int> vec(k, 0);
                    vec[w] = cw;
                    for (auto [y, cy] : rest) vec[y] = cy;
                    vecs.insert(std::move(vec));
                }
            }
        }
        return {vecs.begin(), vecs.end()};
    }
};

struct Result {
    double loss = INF, lower = INF;
    std::vector<std::vector<double>> u;  // P x k
    std::vector<std::vector<int>> M;     // n x k, the multiplicities the usages were fitted to
};

class Scorer {
public:
    Scorer(py::array_t<int, py::array::c_style | py::array::forcecast> seg,
           py::array_t<int, py::array::c_style | py::array::forcecast> variant,
           py::array_t<int, py::array::c_style | py::array::forcecast> total,
           py::array_t<int, py::array::c_style | py::array::forcecast> cluster,
           std::vector<std::vector<std::vector<Call>>> calls,
           std::vector<std::vector<double>> chat,
           std::vector<std::vector<double>> weight,
           double kappa, double eps, int threads, int cap)
        : calls_(std::move(calls)), chat_(std::move(chat)), weight_(std::move(weight)),
          kappa_(kappa), eps_(eps), threads_(std::max(1, threads)), cap_(cap) {
        n_ = seg.size();
        if (variant.ndim() != 2 || variant.shape(0) != n_ || total.ndim() != 2 || total.shape(0) != n_ || cluster.size() != n_)
            throw std::invalid_argument("variant and total must be (SNVs x samples), cluster one entry per SNV");
        P_ = variant.shape(1);
        if (total.shape(1) != P_ || static_cast<int>(calls_.size()) != P_ || static_cast<int>(chat_.size()) != P_ || static_cast<int>(weight_.size()) != P_)
            throw std::invalid_argument("every input needs one entry per sample");
        S_ = calls_[0].size();
        seg_.assign(seg.data(), seg.data() + n_);
        V_.assign(P_, std::vector<int>(n_));
        D_.assign(P_, std::vector<int>(n_));
        for (int i = 0; i < n_; ++i)
            for (int p = 0; p < P_; ++p) {
                V_[p][i] = variant.data()[i * P_ + p];
                D_[p][i] = total.data()[i * P_ + p];
            }
        int C = 0;
        for (int i = 0; i < n_; ++i) C = std::max(C, cluster.data()[i] + 1);
        groups_.assign(C, std::vector<std::vector<int>>(S_));
        for (int i = 0; i < n_; ++i) groups_[cluster.data()[i]][seg_[i]].push_back(i);
    }

    // Scores trees (each (parent, clusters, cn)) from usages u0 (each P x k).
    py::tuple score(const std::vector<std::tuple<std::vector<int>, std::vector<std::vector<int>>, std::vector<std::vector<std::tuple<int, int, int>>>>>& trees,
                    const std::vector<py::array_t<double, py::array::c_style | py::array::forcecast>>& u0, int rounds,
                    int max_iterations, bool multiplicities) {
        const int T = trees.size();
        if (static_cast<int>(u0.size()) != T) throw std::invalid_argument("one u0 per tree");
        std::vector<Tree> ts(T);
        std::vector<std::vector<std::vector<double>>> us(T);
        for (int t = 0; t < T; ++t) {
            auto& [par, cl, cn] = trees[t];
            ts[t] = Tree{par, cl, cn};
            const int k = par.size();
            if (u0[t].ndim() != 2 || u0[t].shape(0) != P_ || u0[t].shape(1) != k) throw std::invalid_argument("u0 must be (samples x clones)");
            us[t].assign(P_, std::vector<double>(k));
            for (int p = 0; p < P_; ++p)
                for (int v = 0; v < k; ++v) us[t][p][v] = u0[t].data()[p * k + v];
        }
        std::vector<Result> res(T);
        {
            py::gil_scoped_release release;
            std::atomic<int> next{0};
            std::exception_ptr error;
            std::mutex error_mutex;
            auto work = [&]() {
                for (int t; (t = next++) < T;) {
                    try {
                        res[t] = score_one(ts[t], us[t], rounds, max_iterations);
                    } catch (...) {
                        std::lock_guard<std::mutex> lock(error_mutex);
                        if (!error) error = std::current_exception();
                        next = T;
                    }
                }
            };
            std::vector<std::thread> pool;
            const int nt = std::min(threads_, T);
            for (int j = 1; j < nt; ++j) pool.emplace_back(work);
            work();
            for (auto& th : pool) th.join();
            if (error) std::rethrow_exception(error);
        }
        py::array_t<double> loss(T), lower(T);
        py::list u;
        for (int t = 0; t < T; ++t) {
            loss.mutable_data()[t] = res[t].loss;
            lower.mutable_data()[t] = res[t].lower;
            const int k = ts[t].parent.size();
            py::array_t<double> a({P_, k});
            const auto& src = res[t].u.empty() ? us[t] : res[t].u;
            for (int p = 0; p < P_; ++p)
                for (int v = 0; v < k; ++v) a.mutable_data()[p * k + v] = src[p][v];
            u.append(a);
        }
        if (!multiplicities) return py::make_tuple(loss, lower, u);
        py::list Ms;
        for (int t = 0; t < T; ++t) {
            if (res[t].M.empty()) {
                Ms.append(py::none());
                continue;
            }
            const int k = ts[t].parent.size();
            py::array_t<int> a({n_, k});
            for (int i = 0; i < n_; ++i)
                for (int v = 0; v < k; ++v) a.mutable_data()[i * k + v] = res[t].M[i][v];
            Ms.append(a);
        }
        return py::make_tuple(loss, lower, u, Ms);
    }

private:
    int n_, P_, S_;
    std::vector<int> seg_;
    std::vector<std::vector<int>> V_, D_;                     // [p][i]
    std::vector<std::vector<std::vector<int>>> groups_;       // [cluster][segment] -> SNVs
    std::vector<std::vector<std::vector<Call>>> calls_;       // [p][s]
    std::vector<std::vector<double>> chat_, weight_;          // [p][s]
    double kappa_, eps_;
    int threads_, cap_;

    struct Group {
        const std::vector<int>* idx;
        int s;
        std::vector<std::vector<int>> O;
    };

    Result score_one(const Tree& tr, const std::vector<std::vector<double>>& u0, int rounds, int max_iterations) const {
        const int k = tr.parent.size();
        // states, in an order where parents come first
        std::vector<std::vector<int>> ch(k);
        for (int v = 1; v < k; ++v) ch[tr.parent[v]].push_back(v);
        std::vector<int> order{0};
        for (size_t j = 0; j < order.size(); ++j)
            for (int c : ch[order[j]]) order.push_back(c);
        std::vector<std::vector<State>> st(S_, std::vector<State>(k, {1, 1}));
        std::vector<std::vector<char>> event(k, std::vector<char>(S_, 0));
        for (int v = 1; v < k; ++v)
            for (auto [s, x, y] : tr.cn[v]) event[v][s] = 1;
        for (size_t j = 1; j < order.size(); ++j) {
            const int v = order[j], p = tr.parent[v];
            for (int s = 0; s < S_; ++s) st[s][v] = st[s][p];
            for (auto [s, x, y] : tr.cn[v]) st[s][v] = {x, y};
        }
        // multiplicity options of every (cluster, segment)
        std::vector<Group> groups;
        for (int w = 1; w < k; ++w)
            for (int c : tr.clusters[w])
                for (int s = 0; s < S_; ++s) {
                    const auto& idx = groups_[c][s];
                    if (idx.empty()) continue;
                    Options opt{st[s], ch, cap_};
                    auto O = opt.of(tr.parent, w, event[w][s]);
                    if (O.empty()) return Result{};  // infinite loss
                    groups.push_back({&idx, s, std::move(O)});
                }
        // squared error of calls of states that no clone has
        std::vector<double> missing(P_, 0.0);
        for (int s = 0; s < S_; ++s) {
            std::set<State> present(st[s].begin(), st[s].end());
            for (int p = 0; p < P_; ++p)
                for (auto [x, y, m] : calls_[p][s])
                    if (!present.count({x, y})) missing[p] += weight_[p][s] * m * m;
        }

        Result r;
        r.u = u0;
        std::vector<std::vector<int>> M(n_, std::vector<int>(k, 0)), prev;
        std::vector<double> lp(P_), l1p(P_);
        for (int round = 0; round < rounds; ++round) {
            for (const Group& g : groups) {
                const int nO = g.O.size();
                std::vector<double> logp(nO * P_), log1mp(nO * P_);
                for (int o = 0; o < nO; ++o)
                    for (int p = 0; p < P_; ++p) {
                        double y = 0.0;
                        for (int v = 0; v < k; ++v) y += g.O[o][v] * r.u[p][v];
                        const double q = std::min(std::max(y / chat_[p][g.s], 1e-12), 1.0 - 1e-12);
                        logp[o * P_ + p] = std::log(q);
                        log1mp[o * P_ + p] = std::log1p(-q);
                    }
                for (int i : *g.idx) {
                    int best = 0;
                    double best_nll = INF;
                    for (int o = 0; o < nO; ++o) {
                        double nll = 0.0;
                        for (int p = 0; p < P_; ++p) nll -= V_[p][i] * logp[o * P_ + p] + (D_[p][i] - V_[p][i]) * log1mp[o * P_ + p];
                        if (nll < best_nll) { best_nll = nll; best = o; }
                    }
                    M[i] = g.O[best];
                }
            }
            if (round > 0 && M == prev) break;
            double L = 0.0, lb = 0.0;
            std::vector<std::vector<double>> u(P_);
            for (int p = 0; p < P_; ++p) {
                CopyNumberADMM::Instance in{tr.parent, st, M, seg_, V_[p], D_[p], calls_[p], chat_[p], weight_[p]};
                CopyNumberADMM::Options opt;
                opt.kappa = kappa_;
                opt.eps = eps_;
                opt.max_iterations = max_iterations;
                opt.initial_usages = r.u[p];
                auto res = CopyNumberADMM::solve(in, opt);
                L += res.loss + missing[p];
                lb += res.lower_bound + missing[p];
                u[p] = std::move(res.usages);
            }
            r.loss = L;
            r.lower = lb;
            r.u = std::move(u);
            prev = M;
        }
        r.M = std::move(prev);
        return r;
    }
};

}  // namespace

PYBIND11_MODULE(_core, m) {
    m.doc() = "Sapling-CN scoring core: the loss of many clone trees at once";
    py::class_<Scorer>(m, "Scorer")
        .def(py::init<py::array_t<int, py::array::c_style | py::array::forcecast>, py::array_t<int, py::array::c_style | py::array::forcecast>,
                      py::array_t<int, py::array::c_style | py::array::forcecast>, py::array_t<int, py::array::c_style | py::array::forcecast>,
                      std::vector<std::vector<std::vector<Call>>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>,
                      double, double, int, int>(),
             py::arg("segment"), py::arg("variant"), py::arg("total"), py::arg("cluster"), py::arg("calls"),
             py::arg("mean_copy_number"), py::arg("segment_weight"), py::arg("kappa") = -1.0, py::arg("eps") = 1.0,
             py::arg("threads") = 8, py::arg("cap") = 256)
        .def("score", &Scorer::score, py::arg("trees"), py::arg("u0"), py::arg("rounds") = 3, py::arg("max_iterations") = 100000, py::arg("multiplicities") = false,
             "Returns (loss, lower bound, usages) of every tree: two arrays and a list of (samples x clones) arrays.\n"
             "With max_iterations, a solve may stop before its gap reaches eps; its loss is still that of feasible usages.\n"
             "With multiplicities, also a list of (SNVs x clones) arrays: the multiplicities the final usages were fitted to\n"
             "(None for a tree with infinite loss).");
}
