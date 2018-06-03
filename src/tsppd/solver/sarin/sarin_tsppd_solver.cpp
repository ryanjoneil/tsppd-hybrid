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

#include <tsppd/solver/sarin/sarin_tsppd_solver.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Solver;
using namespace std;

SarinTSPPDSolver::SarinTSPPDSolver(
    const TSPPDProblem& problem,
    const map<string, string> options,
    TSPSolutionWriter& writer) :
    SarinTSPSolver(problem, options, writer) {

    initialize_tsppd_constraints();
}

void SarinTSPPDSolver::initialize_tsppd_constraints() {
    // +i < -i for all i
    for (auto p : problem.pickup_indices()) {
        auto d = problem.successor_index(p);
        x[d][p].set(GRB_DoubleAttr_UB, 0);
        x[d][p].set(GRB_DoubleAttr_UB, 0);
        y[p][d].set(GRB_DoubleAttr_LB, 1);
    }

    // +0 directly precedes +i for some i
    GRBLinExpr expr = 0;
    for (auto p : problem.pickup_indices())
        expr += x[start_index][p];
    model.addConstr(expr == 1);

    // -0 is directly preceded by -i for some i
    expr = 0;
    for (auto d : problem.delivery_indices())
        expr += x[d][end_index];
    model.addConstr(expr == 1);
}
