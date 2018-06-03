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

#include <tsppd/solver/focacci/propagator/focacci_tsppd_omc_propagator.h>

using namespace Gecode;
using namespace TSPPD::Data;
using namespace TSPPD::Solver;
using namespace std;

FocacciTSPPDOMCPropagator::FocacciTSPPDOMCPropagator(
    Space& home,
    ViewArray<Int::IntView>& next,
    const unsigned int index,
    const TSPPDProblem& problem) :
    Propagator(home),
    next(next),
    index(index),
    problem(problem),
    start_index(0),
    end_index(problem.successor_index(start_index)) {

    next[index].subscribe(home, *this, Int::PC_INT_VAL);
}

FocacciTSPPDOMCPropagator::FocacciTSPPDOMCPropagator(Space& home, FocacciTSPPDOMCPropagator& p) :
    Propagator(home, p),
    next(p.next),
    index(p.index),
    problem(p.problem),
    start_index(p.start_index),
    end_index(p.end_index) {

    next.update(home, p.next);
}

Propagator* FocacciTSPPDOMCPropagator::copy(Space& home) {
    return new (home) FocacciTSPPDOMCPropagator(home, *this);
}

size_t FocacciTSPPDOMCPropagator::dispose(Space& home) {
    next[index].cancel(home, *this, Int::PC_INT_VAL);
    (void) Propagator::dispose(home);
    return sizeof(*this);
}

PropCost FocacciTSPPDOMCPropagator::cost(const Space& home, const ModEventDelta& med) const {
    return PropCost::unary(PropCost::LO);
}

void FocacciTSPPDOMCPropagator::reschedule(Space& home) {
    next[index].reschedule(home, *this, Int::PC_INT_VAL);
}

ExecStatus FocacciTSPPDOMCPropagator::propagate(Space& home, const ModEventDelta& med) {
    unsigned int from = index;
    unsigned int to = next[index].val();

    // Ignore the (-0 +0) arc.
    if (from == end_index && to == start_index)
        return home.ES_SUBSUMED(*this);

    if (problem.has_successor(from) && problem.has_predecessor(to)) {
        int i_m = problem.successor_index(from);
        int j_p = problem.predecessor_index(to);

        // Ignore (+i -i)
        if ((int) from == j_p)
            return home.ES_SUBSUMED(*this);

        // (+i -j) >> !(-i +j) /\ !(+j -i)
        GECODE_ME_CHECK(next[i_m].nq(home, j_p));
        GECODE_ME_CHECK(next[j_p].nq(home, i_m));

    } else if (problem.has_predecessor(from) && problem.has_successor(to)) {
        int i_p = problem.predecessor_index(from);
        int j_m = problem.successor_index(to);

        // (-i +j) -> !(+i -j) /\ !(-j +i)
        GECODE_ME_CHECK(next[i_p].nq(home, j_m));
        GECODE_ME_CHECK(next[j_m].nq(home, i_p));
    }

    return home.ES_SUBSUMED(*this);
}

ExecStatus FocacciTSPPDOMCPropagator::post(
    Space& home,
    ViewArray<Int::IntView>& next,
    const unsigned int index,
    const TSPPDProblem& problem) {

    (void) new (home) FocacciTSPPDOMCPropagator(home, next, index, problem);
    return ES_OK;
}

void TSPPD::Solver::tsppd_omc(Space& home, IntVarArray& next, const TSPPDProblem& problem) {
    GECODE_POST;

    IntVarArgs next_args(next);
    ViewArray<Int::IntView> next_view(home, next_args);

    for (unsigned int index = 0; index < problem.nodes.size(); ++index)
        if (FocacciTSPPDOMCPropagator::post(home, next_view, index, problem) != ES_OK)
            home.fail();
}
