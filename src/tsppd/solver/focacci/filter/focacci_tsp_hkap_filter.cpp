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

#include <tsppd/ap/primal_dual_ap_solver.h>
#include <tsppd/solver/focacci/filter/focacci_tsp_hkap_filter.h>
#include <tsppd/solver/focacci/filter/one_tree/focacci_tsp_one_tree.h>

using namespace Gecode;
using namespace TSPPD::AP;
using namespace TSPPD::Data;
using namespace TSPPD::Solver;
using namespace std;

FocacciTSPHKAPFilter::FocacciTSPHKAPFilter(
    Home home,
    ViewArray<Int::IntView>& next,
    Int::IntView& primal,
    const TSPPDProblem& problem) :
    FocacciTSPHeldKarpFilter(home, next, primal, problem) { }

FocacciTSPHKAPFilter::FocacciTSPHKAPFilter(Space& home, FocacciTSPHKAPFilter& p) :
    FocacciTSPHeldKarpFilter(home, p) { }


Propagator* FocacciTSPHKAPFilter::copy(Space& home) {
    return new (home) FocacciTSPHKAPFilter(home, *this);
}

ExecStatus FocacciTSPHKAPFilter::propagate(Space& home, const ModEventDelta& med) {
    if (primal.assigned() || next.assigned())
        return home.ES_SUBSUMED(*this);

    if (hk_done)
        return ES_FIX;

    OneTree tree(next, problem);
    auto w = tree.bound();

    PrimalDualAPSolver ap(next.size());

    // Objective filtering.
    GECODE_ME_CHECK(primal.gq(home, (int) ceil(w)));

    // Marginal-cost filtering.
    for (int from = 0; from < (int) next.size(); ++from) {
        for (auto to = 0; to < (int) next.size(); ++to) {
            // This only applies to nonbasic feasible arcs in the MST.
            if (!next[from].in(to)) {
                ap.set_bounds({from, to}, 0, 0);
                continue;
            }

            if (tree.has_edge(from, to)) {
                ap.set_bounds({from, to}, 1, 1);
                ap.set_obj({from, to}, 0);
                continue;
            }

            auto mc = tree.marginal_cost(from, to);
            ap.set_obj({from, to}, mc);

            if (w + mc > primal.max())
                GECODE_ME_CHECK(next[from].nq(home, to));
        }
    }


    // Update dual bound based on AP relaxation.
    ap.solve();
    auto z = ap.get_z();

    GECODE_ME_CHECK(primal.gq(home, ((int) ceil(w)) + z));

    // Reduced cost fixing on arcs.
    for (int from = 0; from < (int) next.size(); ++from) {
        for (auto to = next[from].min(); to <= next[from].max(); ++to) {
            if (!next[from].in(to))
                continue;

            auto rc = ap.get_rc({from, to});
            if (w + z + rc > primal.max())
                GECODE_ME_CHECK(next[from].nq(home, to));
        }
    }

    hk_done = true;
    return ES_FIX;
}

ExecStatus FocacciTSPHKAPFilter::post(
    Home home,
    ViewArray<Int::IntView>& next,
    Int::IntView& primal,
    const TSPPDProblem& problem) {

    if (!primal.assigned() && !next.assigned())
        (void) new (home) FocacciTSPHKAPFilter(home, next, primal, problem);
    return ES_OK;
}

void TSPPD::Solver::tsppd_hkap(
    Home home,
    IntVarArray& next,
    IntVar& primal,
    const TSPPDProblem& problem) {

    GECODE_POST;

    IntVarArgs next_args(next);
    ViewArray<Int::IntView> next_view(home, next_args);

    Int::IntView primal_view(primal);

    GECODE_ES_FAIL(FocacciTSPHKAPFilter::post(home, next_view, primal_view, problem));
}
