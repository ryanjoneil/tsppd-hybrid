/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the tsppd program and library for solving           */
/*  Traveling Salesman Problems with Pickup and Delivery. tsppd requires     */
/*  other commercial and open source software to build. tsppd is decribed    */
/*  in the paper "Exact Methods for Solving Traveling Salesman Problems      */
/*  with Pickup and Delivery in Real Time".                                  */
/*                                                                           */
/*  Copyright (C) 2017 Ryan J. O'Neil <roneil1@gmu.edu>                      */
/*                                                                           */
/*  tsppd is distributed under the terms of the ZIB Academic License.        */
/*  You should have received a copy of the ZIB Academic License along with   */
/*  tsppd. See the file LICENSE. If not, email roneil1@gmu.edu.              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iomanip>
#include <memory>

#include <tsppd/solver/focacci/focacci_tsp_space.h>
#include <tsppd/solver/focacci/brancher/focacci_tsp_closest_neighbor_brancher.h>
#include <tsppd/solver/focacci/brancher/focacci_tsp_regret_brancher.h>
#include <tsppd/solver/focacci/brancher/focacci_tsp_sequential_closest_neighbor_brancher.h>
#include <tsppd/solver/focacci/dual/focacci_closest_neighbor_dual.h>

using namespace Gecode;
using namespace TSPPD::Data;
using namespace TSPPD::Solver;
using namespace std;

FocacciTSPSpace::FocacciTSPSpace(const TSPPDProblem& problem) :
    IntMinimizeSpace(),
    problem(problem),
    next(IntVarArray(*this, problem.nodes.size(), 0, problem.nodes.size() - 1)),
    length(IntVar(*this, 0, Int::Limits::max)),
    dual_bound(IntVar(*this, 0, Int::Limits::max)) { }

FocacciTSPSpace::FocacciTSPSpace(FocacciTSPSpace& s) :
    IntMinimizeSpace(s),
    problem(s.problem),
    next(s.next),
    length(s.length),
    dual_bound(s.dual_bound) {

    next.update(*this, s.next);
    length.update(*this, s.length);
    dual_bound.update(*this, s.dual_bound);
}

Gecode::Space* FocacciTSPSpace::copy() {
    return new FocacciTSPSpace(*this);
}

Gecode::IntVar FocacciTSPSpace::cost() const {
    return length;
}

Gecode::IntVar FocacciTSPSpace::dual() const {
    return dual_bound;
}

// This is the solution inspection method for Gist.
void FocacciTSPSpace::print(ostream& out) const {
    out << "length: " << length << endl;
    out << "dual:   " << dual_bound << endl;
    print_int_array(out, next, "next");
}

void FocacciTSPSpace::constrain(const Space& _best) {
    auto best = dynamic_cast<const IntMinimizeSpace*>(&_best);
    if (best == nullptr)
        throw DynamicCastFailed("FocacciTSPSpace::constrain");

    rel(*this, cost() <= best->cost() - 1);
}

void FocacciTSPSpace::initialize_constraints() {
    rel(*this, next[problem.index("-0")] == problem.index("+0"));
    circuit(*this, build_arc_costs(), next, length);
}

void FocacciTSPSpace::initialize_brancher(const FocacciTSPBrancherType brancher_type) {
    IntVarArgs next_args(next);
    ViewArray<Int::IntView> next_view(*this, next_args);

    if (brancher_type == BRANCHER_CN)
        FocacciTSPClosestNeighborBrancher::post(*this, next_view, problem);
    else if (brancher_type == BRANCHER_REGRET)
        FocacciTSPRegretBrancher::post(*this, next_view, problem);
    else if (brancher_type == BRANCHER_SEQ_CN)
        FocacciTSPSequentialClosestNeighborBrancher::post(*this, next_view, problem);
}

void FocacciTSPSpace::initialize_dual(const FocacciTSPDualType dual_type) {
    if (dual_type == DUAL_NONE)
        return;

    rel(*this, length >= dual_bound);

    if (dual_type == DUAL_CN)
        closest_neighbor_dual(*this, next, dual_bound, problem);
}

vector<string> FocacciTSPSpace::solution() const {
    vector<string> s(problem.nodes.size());

    int from = 0;
    for (size_t i = 0; i < problem.nodes.size(); ++i) {
        int to = next[from].val();
        s[i] = problem.nodes[from];
        from = to;
    }

    return s;
}

IntArgs FocacciTSPSpace::build_arc_costs() const {
    IntArgs arc_costs(problem.nodes.size() * problem.nodes.size());
    for (size_t from = 0; from < problem.nodes.size(); ++from) {
        for (size_t to = 0; to < problem.nodes.size(); ++to) {
            auto index = (from * problem.nodes.size()) + to;

            if (from == to)
                arc_costs[index] = 0;
            else
                arc_costs[index] = problem.cost(problem.nodes[from], problem.nodes[to]);
        }
    }
    return arc_costs;
}

void FocacciTSPSpace::print_int_array(ostream& out, IntVarArray array, string label) const {
    out << label << ": ";
    for (int index = 0; index < array.size(); ++index) {
        if (index > 0)
            out << setfill(' ') << setw(label.length() + 3);

        auto node = problem.nodes[index];
        auto var = array[index];
        out << "[" << node << "] ";
        for (auto i = var.min(); i <= var.max(); ++i)
            if (var.in(i))
                out << problem.nodes[i] << " ";

        out << endl;
    }
}
