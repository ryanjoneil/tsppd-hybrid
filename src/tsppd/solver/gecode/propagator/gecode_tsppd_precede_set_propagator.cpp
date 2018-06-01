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

#include <tsppd/solver/gecode/propagator/gecode_tsppd_precede_set_propagator.h>

using namespace Gecode;
using namespace TSPPD::Data;
using namespace TSPPD::Solver;
using namespace std;

GecodeTSPPDPrecedeSetPropagator::GecodeTSPPDPrecedeSetPropagator(
    Home home,
    ViewArray<Int::IntView>& next,
    ViewArray<Set::SetView>& pred,
    ViewArray<Set::SetView>& succ,
    const unsigned int index,
    const TSPPDProblem& problem) :
    Propagator(home),
    next(next),
    pred(pred),
    succ(succ),
    index(index),
    problem(problem),
    start_index(0),
    end_index(problem.successor_index(start_index)) {

    next[index].subscribe(home, *this, Int::PC_INT_VAL);
}

GecodeTSPPDPrecedeSetPropagator::GecodeTSPPDPrecedeSetPropagator(Space& home, GecodeTSPPDPrecedeSetPropagator& p) :
    Propagator(home, p),
    next(p.next),
    pred(p.pred),
    succ(p.succ),
    index(p.index),
    problem(p.problem),
    start_index(p.start_index),
    end_index(p.end_index) {

    next.update(home, p.next);
    pred.update(home, p.pred);
    succ.update(home, p.succ);
}

Propagator* GecodeTSPPDPrecedeSetPropagator::copy(Space& home) {
    return new (home) GecodeTSPPDPrecedeSetPropagator(home, *this);
}

size_t GecodeTSPPDPrecedeSetPropagator::dispose(Space& home) {
    next[index].cancel(home, *this, Int::PC_INT_VAL);
    (void) Propagator::dispose(home);
    return sizeof(*this);
}

PropCost GecodeTSPPDPrecedeSetPropagator::cost(const Space& home, const ModEventDelta& med) const {
    return PropCost::crazy(PropCost::LO, (unsigned int) next.size());
}

void GecodeTSPPDPrecedeSetPropagator::reschedule(Space& home) {
    next[index].reschedule(home, *this, Int::PC_INT_VAL);
}

ExecStatus GecodeTSPPDPrecedeSetPropagator::propagate(Space& home, const ModEventDelta& med) {
    unsigned int from = index;
    unsigned int to = next[index].val();

    // Ignore the (-0 +0) arc.
    if (from == end_index && to == start_index)
        return home.ES_SUBSUMED(*this);

    // succ[from] = succ[to] + {to}
    rel(home, succ[to], SOT_UNION, IntSet(to, to), SRT_EQ, succ[from]);

    // pred[to] = pred[from] + [from]
    rel(home, pred[from], SOT_UNION, IntSet(from, from), SRT_EQ, pred[to]);

    GECODE_ME_CHECK(precede(home, from, to));

    return home.ES_SUBSUMED(*this);
}

ExecStatus GecodeTSPPDPrecedeSetPropagator::post(
    Home home,
    ViewArray<Int::IntView>& next,
    ViewArray<Set::SetView>& pred,
    ViewArray<Set::SetView>& succ,
    const unsigned int index,
    const TSPPDProblem& problem) {

    (void) new (home) GecodeTSPPDPrecedeSetPropagator(home, next, pred, succ, index, problem);
    return ES_OK;
}

ExecStatus GecodeTSPPDPrecedeSetPropagator::precede(Space& home, const unsigned int i, const unsigned int j) {
    GECODE_ME_CHECK(before(home, i, j));

    if (problem.has_successor(i)) {
        auto i_p = i;
        auto i_m = problem.successor_index(i);

        if (problem.has_successor(j)) {
            // (+i +j)
            auto j_p = j;
            auto j_m = problem.successor_index(j);

            // -i after +j
            GECODE_ME_CHECK(after(home, i_m, j_p));

            // -j after +i
            GECODE_ME_CHECK(after(home, j_m, i_p));

        } else if (problem.has_predecessor(j)) {
            // (+i -j)
            auto j_p = problem.predecessor_index(j);
            auto j_m = j;

            // +j before {+i, -i}
            GECODE_ME_CHECK(before(home, j_p, i_p));
            GECODE_ME_CHECK(before(home, j_p, i_m));

            // -i after -j
            GECODE_ME_CHECK(after(home, i_m, j_m));

            if (i_p != j_p) {
                // !( -i +j )
                GECODE_ME_CHECK(next[i_m].nq(home, (int) j_p));

                // !( +j -i )
                GECODE_ME_CHECK(next[j_p].nq(home, (int) i_m));
            }
        }

    } else if (problem.has_predecessor(i)) {
        auto i_p = problem.predecessor_index(i);
        auto i_m = i;

        if (problem.has_successor(j)) {
            // (-i +j)
            auto j_p = j;
            auto j_m = problem.successor_index(j);

            // +i before {+j, -j}
            GECODE_ME_CHECK(before(home, i_p, j_p));
            GECODE_ME_CHECK(before(home, i_p, j_m));

            // -j after -i
            GECODE_ME_CHECK(after(home, j_m, i_m));

            if (i_p != j_p) {
                // !( -j +i )
                GECODE_ME_CHECK(next[j_m].nq(home, (int) i_p));

                // !( +i -j )
                GECODE_ME_CHECK(next[i_p].nq(home, (int) j_m));
            }

        } else if (problem.has_predecessor(j)) {
            // (-i -j)
            auto j_p = problem.predecessor_index(j);
            auto j_m = j;

            // +i before -j
            GECODE_ME_CHECK(before(home, i_p, j_m));

            // +j before -i
            GECODE_ME_CHECK(before(home, j_p, i_m));
        }
    }

    return ES_OK;
}

ExecStatus GecodeTSPPDPrecedeSetPropagator::before(Space& home, const unsigned int i, const unsigned int j) {
    if (i == j)
        return ES_OK;

    GECODE_ME_CHECK(next[j].nq(home, (int) i));

    // i before j
    if (i != end_index)
        GECODE_ME_CHECK(succ[i].include(home, j));

    if (i != start_index)
        GECODE_ME_CHECK(pred[i].exclude(home, j));

    // j after i
    if (j != start_index)
        GECODE_ME_CHECK(pred[j].include(home, i));

    if (j != end_index)
        GECODE_ME_CHECK(succ[j].exclude(home, i));

    return ES_OK;
}

ExecStatus GecodeTSPPDPrecedeSetPropagator::after(Space& home, const unsigned int i, const unsigned int j) {
    // i after j = j before i
    return before(home, j, i);
}

void TSPPD::Solver::tsppd_precede_set(Home home, IntVarArray& next, const TSPPDProblem& problem) {
    GECODE_POST;

    // Construct predecessor and successor sets for each node.
    SetVarArray pred(home, problem.nodes.size(), IntSet::empty, IntSet(0, problem.nodes.size() - 1));
    SetVarArray succ(home, problem.nodes.size(), IntSet::empty, IntSet(0, problem.nodes.size() - 1));

    IntVarArgs next_args(next);
    ViewArray<Int::IntView> next_view(home, next_args);

    SetVarArgs pred_args(pred);
    ViewArray<Set::SetView> pred_view(home, pred_args);

    SetVarArgs succ_args(succ);
    ViewArray<Set::SetView> succ_view(home, succ_args);

    unsigned int end_index = problem.successor_index(0);
    for (unsigned int index = 0; index < problem.nodes.size(); ++index) {
        if (index == 0)
            // The start node has nothing before it and all other nodes after it.
            rel(home, pred[index] == IntSet::empty);

        else if (index == end_index)
            // The End node has nothing after it and all other nodes before it.
            rel(home, succ[index] == IntSet::empty);

        if (problem.has_successor(index)) {
            // Pickups cannot be preceded and must be succeeded by their deliveries.
            auto delivery = problem.successor_index(index);
            dom(home, pred[index], SRT_DISJ, delivery);
            dom(home, succ[index], SRT_SUP, delivery);
            rel(home, succ[index], SRT_SUP, succ[delivery]);

        } else if (problem.has_predecessor(index)) {
            // Deliveries cannot be succeeded and must be preceded by their pickups.
            auto pickup = problem.predecessor_index(index);
            dom(home, pred[index], SRT_SUP, pickup);
            dom(home, succ[index], SRT_DISJ, pickup);
            rel(home, pred[index], SRT_SUP, pred[pickup]);
        }

        GECODE_ES_FAIL(GecodeTSPPDPrecedeSetPropagator::post(home, next_view, pred_view, succ_view, index, problem));
    }
}