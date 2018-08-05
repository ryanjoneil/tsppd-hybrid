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

#include <tsppd/solver/focacci/filter/focacci_tsp_assignment_filter.h>

using namespace Gecode;
using namespace TSPPD::AP;
using namespace TSPPD::Data;
using namespace TSPPD::Solver;
using namespace std;

FocacciTSPAssignmentFilter::FocacciTSPAssignmentFilter(
    Home home,
    ViewArray<Int::IntView>& next,
    Int::IntView& primal,
    const TSPPDProblem& problem) :
    Propagator(home),
    next(next),
    primal(primal),
    problem(problem),
    ap(PrimalDualAPSolver(next.size())),
    unassigned() {

    // Set arc costs on AP relaxation.
    for (int from = 0; from < next.size(); ++from) {
        for (int to = 0; to < next.size(); ++to) {
            ap.set_obj({from, to}, problem.cost(from, to));

            if (!next[from].in(to))
                ap.set_bounds({from, to}, 0, 0);

            else if (next[from].assigned() && next[from].val() == to)
                ap.set_bounds({from, to}, 1, 1);

            else
                unassigned.push_back({from, to});
        }
    }

    next.subscribe(home, *this, Int::PC_INT_DOM);
    home.notice(*this, AP_DISPOSE);
}

FocacciTSPAssignmentFilter::FocacciTSPAssignmentFilter(Space& home, FocacciTSPAssignmentFilter& p) :
    Propagator(home, p),
    next(p.next),
    primal(p.primal),
    problem(p.problem),
    ap(p.ap),
    unassigned(p.unassigned) {

    next.update(home, p.next);
    primal.update(home, p.primal);
}

Propagator* FocacciTSPAssignmentFilter::copy(Space& home) {
    return new (home) FocacciTSPAssignmentFilter(home, *this);
}

size_t FocacciTSPAssignmentFilter::dispose(Space& home) {
    home.ignore(*this, AP_DISPOSE);
    next.cancel(home, *this, Int::PC_INT_DOM);
    ap.~PrimalDualAPSolver();
    unassigned.~vector<pair<int, int>>();
    (void) Propagator::dispose(home);
    return sizeof(*this);
}

PropCost FocacciTSPAssignmentFilter::cost(const Space& home, const ModEventDelta& med) const {
    return PropCost::quadratic(PropCost::HI, next.size());
}

void FocacciTSPAssignmentFilter::reschedule(Space& home) {
    next.reschedule(home, *this, Int::PC_INT_DOM);
}

ExecStatus FocacciTSPAssignmentFilter::propagate(Space& home, const ModEventDelta& med) {
    if (primal.assigned() || next.assigned() || unassigned.empty())
        return home.ES_SUBSUMED(*this);

    // Update variable bounds.
    vector<pair<int, int>> new_unassigned;

    for (auto indices : unassigned) {
        auto from = indices.first;
        auto to = indices.second;

        bool lb = false;
        bool ub = true;

        if (next[from].assigned()) {
            lb = (next[from].val() == to);
            ub = lb;
        } else if (!next[from].in(to))
            ub = false;
        else
            new_unassigned.push_back({from, to});

        ap.set_bounds({from, to}, lb, ub);
    }

    if (new_unassigned.size() > 0)
        unassigned = new_unassigned;
    else
        return home.ES_SUBSUMED(*this);

    // Make sure the solution succeeds.
    if (!ap.solve())
        return ES_FAILED;

    // Update dual bound based on AP relaxation.
    auto z = ap.get_z();
    GECODE_ME_CHECK(primal.gq(home, z));

    // Reduced cost fixing on arcs.
    for (auto indices : unassigned) {
        auto from = indices.first;
        auto to = indices.second;

        auto c = ap.get_rc({from, to});
        if (z + c > primal.max())
            GECODE_ME_CHECK(next[from].nq(home, to));
    }

    return ES_FIX;
}

ExecStatus FocacciTSPAssignmentFilter::post(
    Home home,
    ViewArray<Int::IntView>& next,
    Int::IntView& primal,
    const TSPPDProblem& problem) {

    if (!primal.assigned() && !next.assigned())
        (void) new (home) FocacciTSPAssignmentFilter(home, next, primal, problem);
    return ES_OK;
}

void TSPPD::Solver::tsppd_assignment(
    Home home,
    IntVarArray& next,
    IntVar& primal,
    const TSPPDProblem& problem) {

    GECODE_POST;

    IntVarArgs next_args(next);
    ViewArray<Int::IntView> next_view(home, next_args);

    Int::IntView primal_view(primal);

    GECODE_ES_FAIL(FocacciTSPAssignmentFilter::post(home, next_view, primal_view, problem));
}
