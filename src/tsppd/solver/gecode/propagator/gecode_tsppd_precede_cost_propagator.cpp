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

#include <tsppd/solver/gecode/propagator/gecode_tsppd_precede_cost_propagator.h>

using namespace Gecode;
using namespace TSPPD::Data;
using namespace TSPPD::Solver;
using namespace std;

GecodeTSPPDPrecedeCostPropagator::GecodeTSPPDPrecedeCostPropagator(
    Home home,
    ViewArray<Int::IntView>& next,
    ViewArray<Int::IntView>& node_cost,
    const unsigned int index,
    const TSPPDProblem& problem) :
    Propagator(home),
    next(next),
    node_cost(node_cost),
    index(index),
    problem(problem),
    start_index(0),
    end_index(problem.successor_index(start_index)) {

    next[index].subscribe(home, *this, Int::PC_INT_VAL);
    node_cost[index].subscribe(home, *this, Int::PC_INT_DOM);
}

GecodeTSPPDPrecedeCostPropagator::GecodeTSPPDPrecedeCostPropagator(Space& home, GecodeTSPPDPrecedeCostPropagator& p) :
    Propagator(home, p),
    next(p.next),
    node_cost(p.node_cost),
    index(p.index),
    problem(p.problem),
    start_index(p.start_index),
    end_index(p.end_index) {

    next.update(home, p.next);
    node_cost.update(home, p.node_cost);
}

Propagator* GecodeTSPPDPrecedeCostPropagator::copy(Space& home) {
    return new (home) GecodeTSPPDPrecedeCostPropagator(home, *this);
}

size_t GecodeTSPPDPrecedeCostPropagator::dispose(Space& home) {
    next[index].cancel(home, *this, Int::PC_INT_VAL);
    node_cost[index].cancel(home, *this, Int::PC_INT_DOM);

    (void) Propagator::dispose(home);
    return sizeof(*this);
}

PropCost GecodeTSPPDPrecedeCostPropagator::cost(const Space& home, const ModEventDelta& med) const {
    return PropCost::unary(PropCost::LO);
}

void GecodeTSPPDPrecedeCostPropagator::reschedule(Space& home) {
    node_cost[index].reschedule(home, *this, Int::PC_INT_DOM);
    next[index].reschedule(home, *this, Int::PC_INT_VAL);
}

ExecStatus GecodeTSPPDPrecedeCostPropagator::propagate(Space& home, const ModEventDelta& med) {
    if (!next[index].assigned())
        return ES_FIX;

    unsigned int from = index;
    unsigned int to = next[index].val();

    // Ignore the (-0 +0) arc.
    if (from == end_index && to == start_index)
        return home.ES_SUBSUMED(*this);

    auto c = problem.cost(from, to);
    if (node_cost[from].assigned())
        GECODE_ME_CHECK(node_cost[to].eq(home, node_cost[from].val() + c));
    else if (node_cost[to].assigned())
        GECODE_ME_CHECK(node_cost[from].eq(home, node_cost[to].val() - c));
    else {
        GECODE_ME_CHECK(node_cost[to].gq(home, node_cost[from].min() + c));
        GECODE_ME_CHECK(node_cost[to].lq(home, node_cost[from].max() + c));

        GECODE_ME_CHECK(node_cost[from].gq(home, node_cost[to].min() - c));
        GECODE_ME_CHECK(node_cost[from].lq(home, node_cost[to].max() - c));
    }

    if (node_cost[from].assigned() && node_cost[to].assigned())
        return home.ES_SUBSUMED(*this);

    return ES_FIX;
}

ExecStatus GecodeTSPPDPrecedeCostPropagator::post(
    Home home,
    ViewArray<Int::IntView>& next,
    ViewArray<Int::IntView>& node_cost,
    const unsigned int index,
    const TSPPDProblem& problem) {

    (void) new (home) GecodeTSPPDPrecedeCostPropagator(home, next, node_cost, index, problem);
    return ES_OK;
}

void TSPPD::Solver::tsppd_precede_cost(
    Home home,
    IntVar& length,
    IntVarArray& next,
    const TSPPDProblem& problem) {

    GECODE_POST;

    IntVarArgs next_args(next);
    ViewArray<Int::IntView> next_view(home, next_args);

    IntVarArray node_cost(home, problem.nodes.size(), 0, Int::Limits::max);
    IntVarArgs node_cost_args(node_cost);
    ViewArray<Int::IntView> node_cost_view(home, node_cost_args);

    unsigned int end_index = problem.successor_index(0);
    for (unsigned int index = 0; index < problem.nodes.size(); ++index) {
        rel(home, length >= node_cost[index]);

        if (index == 0)
            rel(home, node_cost[index] == 0);
        else if (index == end_index)
            continue;
        else if (problem.has_predecessor(index)) {
            auto pickup = problem.predecessor_index(index);
            rel(home, node_cost[index] >= node_cost[pickup] + problem.cost(pickup, index));
        }

        GECODE_ES_FAIL(GecodeTSPPDPrecedeCostPropagator::post(home, next_view, node_cost_view, index, problem));
    }
}