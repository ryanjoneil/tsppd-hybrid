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

#include <limits>

#include <tsppd/solver/focacci/brancher/focacci_tsp_closest_neighbor_brancher.h>

using namespace Gecode;
using namespace TSPPD::Data;
using namespace TSPPD::Solver;
using namespace std;

FocacciTSPClosestNeighborBrancher::FocacciTSPClosestNeighborBrancher(
    Home home,
    ViewArray<Int::IntView>& next,
    const TSPPDProblem& problem) :
    FocacciTSPBrancher(home, next, problem) {

    home.notice(*this, AP_DISPOSE);
}

FocacciTSPClosestNeighborBrancher::FocacciTSPClosestNeighborBrancher(
    Space& home,
    FocacciTSPClosestNeighborBrancher& b) :
    FocacciTSPBrancher(home, b) { }

Actor* FocacciTSPClosestNeighborBrancher::copy(Space& home) {
    return new (home) FocacciTSPClosestNeighborBrancher(home, *this);
}

size_t FocacciTSPClosestNeighborBrancher::dispose(Gecode::Space& home) {
    home.ignore(*this, AP_DISPOSE);
    (void) FocacciTSPBrancher::dispose(home);
    return sizeof(*this);
}

Choice* FocacciTSPClosestNeighborBrancher::choice(Space& home) {    // Scan for max regret
    TSPPDArc best_arc;

    for (int from = 0; from < next.size(); ++from) {
        if (next[from].assigned() || next[from].size() < 2)
            continue;

        auto arc_index = closest_feasible_arc_index(from, indexes[from]);
        auto arc = problem.arc(from, arc_index);

        indexes[from] = arc_index;

        if (arc.cost < best_arc.cost)
            best_arc = arc;
    }

    return new FocacciTSPBranchChoice(*this, best_arc.from_index, best_arc.to_index);
}

void FocacciTSPClosestNeighborBrancher::post(
    Home home,
    ViewArray<Int::IntView>& next,
    const TSPPDProblem& problem) {

    (void) new (home) FocacciTSPClosestNeighborBrancher(home, next, problem);
}
