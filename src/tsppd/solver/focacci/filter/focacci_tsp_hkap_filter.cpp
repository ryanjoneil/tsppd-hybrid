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

#include <tsppd/solver/focacci/filter/focacci_tsp_hkap_filter.h>

using namespace Gecode;
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
    auto status = FocacciTSPHeldKarpFilter::propagate(home, med);
    if (status != ES_FIX)
        return status;

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
