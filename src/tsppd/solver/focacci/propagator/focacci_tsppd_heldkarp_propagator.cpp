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

#include <cmath>

#include <tsppd/solver/focacci/one_tree/focacci_tsppd_one_tree.h>
#include <tsppd/solver/focacci/propagator/focacci_tsppd_heldkarp_propagator.h>

using namespace Gecode;
using namespace TSPPD::AP;
using namespace TSPPD::Data;
using namespace TSPPD::Solver;
using namespace std;

FocacciTSPPDHeldKarpPropagator::FocacciTSPPDHeldKarpPropagator(
    Home home,
    ViewArray<Int::IntView>& next,
    Int::IntView& primal,
    const TSPPDProblem& problem) :
    Propagator(home),
    next(next),
    primal(primal),
    problem(problem),
    start_index(problem.index("+0")),
    end_index(problem.index("-0")) {

    next.subscribe(home, *this, Int::PC_INT_DOM);
    home.notice(*this, AP_DISPOSE);
}

FocacciTSPPDHeldKarpPropagator::FocacciTSPPDHeldKarpPropagator(Space& home, FocacciTSPPDHeldKarpPropagator& p) :
    Propagator(home, p),
    next(p.next),
    primal(p.primal),
    problem(p.problem),
    start_index(p.start_index),
    end_index(p.end_index) {

    next.update(home, p.next);
    primal.update(home, p.primal);
}

Propagator* FocacciTSPPDHeldKarpPropagator::copy(Space& home) {
    return new (home) FocacciTSPPDHeldKarpPropagator(home, *this);
}

size_t FocacciTSPPDHeldKarpPropagator::dispose(Space& home) {
    home.ignore(*this, AP_DISPOSE);
    next.cancel(home, *this, Int::PC_INT_DOM);
    (void) Propagator::dispose(home);
    return sizeof(*this);
}

PropCost FocacciTSPPDHeldKarpPropagator::cost(const Space& home, const ModEventDelta& med) const {
    return PropCost::crazy(PropCost::HI, next.size());
}

void FocacciTSPPDHeldKarpPropagator::reschedule(Space& home) {
    next.reschedule(home, *this, Int::PC_INT_DOM);
}

ExecStatus FocacciTSPPDHeldKarpPropagator::propagate(Space& home, const ModEventDelta& med) {
    if (primal.assigned() || next.assigned())
        return home.ES_SUBSUMED(*this);

    OneTree tree(next, problem);
    auto w = tree.bound();

    // Objective filtering.
    GECODE_ME_CHECK(primal.gq(home, (int) ceil(w)));

    // Marginal-cost filtering.
    for (int from = 0; from < (int) next.size(); ++from) {
        for (auto to = next[from].min(); to <= next[from].max(); ++to) {
            // This only applies to nonbasic feasible arcs in the MST.
            if (!(next[from].in(to) && !tree.has_edge(from, to)))
                continue;

            if (w + tree.marginal_cost(from, to) > primal.max())
                GECODE_ME_CHECK(next[from].nq(home, to));
        }
    }

    return ES_FIX;
}

ExecStatus FocacciTSPPDHeldKarpPropagator::post(
    Home home,
    ViewArray<Int::IntView>& next,
    Int::IntView& primal,
    const TSPPDProblem& problem) {

    if (!primal.assigned() && !next.assigned())
        (void) new (home) FocacciTSPPDHeldKarpPropagator(home, next, primal, problem);
    return ES_OK;
}

void TSPPD::Solver::tsppd_heldkarp(
    Home home,
    IntVarArray& next,
    IntVar& primal,
    const TSPPDProblem& problem) {

    GECODE_POST;

    IntVarArgs next_args(next);
    ViewArray<Int::IntView> next_view(home, next_args);

    Int::IntView primal_view(primal);

    GECODE_ES_FAIL(FocacciTSPPDHeldKarpPropagator::post(home, next_view, primal_view, problem));
}
