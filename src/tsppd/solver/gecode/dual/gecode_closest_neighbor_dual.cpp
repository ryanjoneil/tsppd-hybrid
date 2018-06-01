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

#include <tsppd/solver/gecode/dual/gecode_closest_neighbor_dual.h>

using namespace Gecode;
using namespace TSPPD::Data;
using namespace TSPPD::Solver;
using namespace std;

GecodeClosestNeighborDual::GecodeClosestNeighborDual(
    Space& home,
    Int::IntView next,
    Int::IntView closest_cost,
    const unsigned int node_index,
    const TSPPDProblem& problem) :
    Propagator(home),
    next(next),
    closest_cost(closest_cost),
    node_index(node_index),
    problem(problem),
    arc_index(0) {

    next.subscribe(home, *this, Int::PC_INT_DOM);
}

GecodeClosestNeighborDual::GecodeClosestNeighborDual(
    Space& home,
    GecodeClosestNeighborDual& p) :
    Propagator(home, p),
    next(p.next),
    closest_cost(p.closest_cost),
    node_index(p.node_index),
    problem(p.problem),
    arc_index(p.arc_index) {

    next.update(home, p.next);
    closest_cost.update(home, p.closest_cost);
}

Propagator* GecodeClosestNeighborDual::copy(Space& home) {
    return new (home) GecodeClosestNeighborDual(home, *this);
}

size_t GecodeClosestNeighborDual::dispose(Space& home) {
    next.cancel(home, *this, Int::PC_INT_DOM);
    (void) Propagator::dispose(home);
    return sizeof(*this);
}

PropCost GecodeClosestNeighborDual::cost(const Space& home, const ModEventDelta& med) const {
    return PropCost::linear(PropCost::HI, ((int) problem.arcs_size(node_index)) - arc_index - 1);
}

void GecodeClosestNeighborDual::reschedule(Space& home) {
    next.reschedule(home, *this, Int::PC_INT_DOM);
}

ExecStatus GecodeClosestNeighborDual::propagate(Space& home, const ModEventDelta& med) {
    while (arc_index < problem.arcs_size(node_index)) {
        auto arc = problem.arc(node_index, arc_index);
        if (next.in((int) arc.to_index)) {
            if (closest_cost.gq(home, arc.cost) == Int::ME_INT_FAILED)
                return ES_FAILED;
            return ES_NOFIX;
        }
        ++arc_index;
    }

    return home.ES_SUBSUMED(*this);
}

ExecStatus GecodeClosestNeighborDual::post(
    Space& home,
    Int::IntView next,
    Int::IntView closest_cost,
    const unsigned int node_index,
    const TSPPDProblem& problem) {

    (void) new (home) GecodeClosestNeighborDual(home, next, closest_cost, node_index, problem);
    return ES_OK;
}

void TSPPD::Solver::closest_neighbor_dual(
    Space& home,
    IntVarArray next,
    IntVar dual,
    const TSPPDProblem& problem) {

    GECODE_POST;

    IntVarArgs costs;
    for (unsigned int node_index = 0; node_index < problem.nodes.size(); ++node_index) {
        // Closest successor cost
        Gecode::IntVar closest_next_cost(home, 0,  Int::Limits::max);
        costs << closest_next_cost;

        Int::IntView next_view(next[node_index]);
        Int::IntView closest_next_cost_view(closest_next_cost);
        auto result = GecodeClosestNeighborDual::post(
            home,
            next_view,
            closest_next_cost_view,
            node_index,
            problem
        );
        if (result != ES_OK)
            home.fail();
    }

    rel(home, dual == sum(costs));
}
