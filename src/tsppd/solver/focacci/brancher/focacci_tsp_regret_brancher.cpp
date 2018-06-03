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

#include <tsppd/solver/focacci/brancher/focacci_tsp_regret_brancher.h>

using namespace Gecode;
using namespace TSPPD::Data;
using namespace TSPPD::Solver;
using namespace std;

FocacciTSPRegretBrancher::FocacciTSPRegretBrancher(
    Home home,
    ViewArray<Int::IntView>& next,
    const TSPPDProblem& problem) :
    FocacciTSPBrancher(home, next, problem) {

    home.notice(*this, AP_DISPOSE);
}

FocacciTSPRegretBrancher::FocacciTSPRegretBrancher(Space& home, FocacciTSPRegretBrancher& b) :
    FocacciTSPBrancher(home, b) { }

Actor* FocacciTSPRegretBrancher::copy(Space& home) {
    return new (home) FocacciTSPRegretBrancher(home, *this);
}

size_t FocacciTSPRegretBrancher::dispose(Gecode::Space& home) {
    home.ignore(*this, AP_DISPOSE);
    (void) FocacciTSPBrancher::dispose(home);
    return sizeof(*this);
}

Choice* FocacciTSPRegretBrancher::choice(Space& home) {
    // Scan for max regret
    int max_regret = 0;
    int max_regret_from = 0;
    int max_regret_to = 0;

    for (int from = 0; from < next.size(); ++from) {
        if (next[from].assigned() || next[from].size() < 2)
            continue;

        auto index_1 = closest_feasible_arc_index(from, indexes[from]);
        indexes[from] = index_1;
        if (index_1 >= problem.arcs_size(from))
            continue;

        auto index_2 = closest_feasible_arc_index(from, index_1 + 1);
        if (index_2 >= problem.arcs_size(from)) {
            indexes[from] = problem.arcs_size(from);
            continue;
        }

        auto arc_1 = problem.arc(from, index_1);
        auto arc_2 = problem.arc(from, index_2);

        auto regret = arc_2.cost - arc_1.cost;
        if (max_regret == 0 || regret > max_regret) {
            max_regret = regret;
            max_regret_from = from;
            max_regret_to = arc_1.to_index;;
        }
    }

    return new FocacciTSPBranchChoice(*this, max_regret_from, max_regret_to);
}

void FocacciTSPRegretBrancher::post(
    Home home,
    ViewArray<Int::IntView>& next,
    const TSPPDProblem& problem) {

    (void) new (home) FocacciTSPRegretBrancher(home, next, problem);
}
