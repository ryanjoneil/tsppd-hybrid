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

#include <tsppd/solver/ap/ap_atsppd_callback.h>
#include <tsppd/solver/ap/ap_atsppd_solver.h>
#include <tsppd/util/exception.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Solver;
using namespace TSPPD::Util;
using namespace std;

APATSPPDSolver::APATSPPDSolver(
    const TSPPDProblem& problem,
    const map<string, string> options,
    TSPSolutionWriter& writer) :
    APATSPSolver(problem, options, writer) {

    initialize_atsppd_variables();
}

TSPPDSolution APATSPPDSolver::solve() {
    configure_solver();

    if (sec == ATSP_SEC_OTHER)
        throw TSPPDException("sec can be either subtour or cutset");

    APATSPPDCallback callback(problem, x, sec, writer);
    model.setCallback(&callback);

    model.optimize();
    return solution();
}

void APATSPPDSolver::initialize_atsppd_variables() {
    for (auto p : problem.pickup_indices()) {
        // ~(+0 -i)
        auto d = problem.successor_index(p);
        x[start_index][d].set(GRB_DoubleAttr_UB, 0);

        // ~(+i -0)
        x[p][end_index].set(GRB_DoubleAttr_UB, 0);

        // ~(-i +i)
        x[d][p].set(GRB_DoubleAttr_UB, 0);
    }
}