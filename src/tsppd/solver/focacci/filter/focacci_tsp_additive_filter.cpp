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

#include <tsppd/solver/focacci/filter/focacci_tsp_additive_filter.h>

using namespace Gecode;
using namespace TSPPD::AP;
using namespace TSPPD::Data;
using namespace TSPPD::Solver;
using namespace std;

FocacciTSPAdditiveFilter::FocacciTSPAdditiveFilter(
    Home home,
    ViewArray<Int::IntView>& next,
    Int::IntView& primal,
    const TSPPDProblem& problem) :
    FocacciTSPAssignmentFilter(home, next, primal, problem) { }

FocacciTSPAdditiveFilter::FocacciTSPAdditiveFilter(Space& home, FocacciTSPAdditiveFilter& p) :
    FocacciTSPAssignmentFilter(home, p) { }


Propagator* FocacciTSPAdditiveFilter::copy(Space& home) {
    return new (home) FocacciTSPAdditiveFilter(home, *this);
}

ExecStatus FocacciTSPAdditiveFilter::propagate(Space& home, const ModEventDelta& med) {
    auto status = FocacciTSPAssignmentFilter::propagate(home, med);
    if (status != ES_FIX)
        return status;

    if (hk_done)
        return ES_FIX;

    OneTree tree(next, problem, &ap);
    auto z = ap.get_z();
    auto w = tree.bound();

    // Objective filtering.
    GECODE_ME_CHECK(primal.gq(home, z + (int) ceil(w)));

    // Marginal-cost filtering.
    for (int from = 0; from < (int) next.size(); ++from) {
        for (auto to = next[from].min(); to <= next[from].max(); ++to) {
            // This only applies to nonbasic feasible arcs in the MST.
            if (!next[from].in(to) || tree.has_edge(from, to))
                continue;

            auto rc = ap.get_rc({from, to});
            auto mc = tree.marginal_cost(from, to);
            if (z + w + rc + mc > primal.max())
                GECODE_ME_CHECK(next[from].nq(home, to));
        }
    }

    hk_done = true;
    return ES_FIX;
}

ExecStatus FocacciTSPAdditiveFilter::post(
    Home home,
    ViewArray<Int::IntView>& next,
    Int::IntView& primal,
    const TSPPDProblem& problem) {

    if (!primal.assigned() && !next.assigned())
        (void) new (home) FocacciTSPAdditiveFilter(home, next, primal, problem);
    return ES_OK;
}

void TSPPD::Solver::tsppd_additive(
    Home home,
    IntVarArray& next,
    IntVar& primal,
    const TSPPDProblem& problem) {

    GECODE_POST;

    IntVarArgs next_args(next);
    ViewArray<Int::IntView> next_view(home, next_args);

    Int::IntView primal_view(primal);

    GECODE_ES_FAIL(FocacciTSPAdditiveFilter::post(home, next_view, primal_view, problem));
}
