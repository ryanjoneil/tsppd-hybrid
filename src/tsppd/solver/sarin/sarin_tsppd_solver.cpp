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
#include <tsppd/util/exception.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Solver;
using namespace TSPPD::Util;
using namespace std;

SarinTSPPDSolver::SarinTSPPDSolver(
    const TSPPDProblem& problem,
    const map<string, string> options,
    TSPSolutionWriter& writer) :
    SarinTSPSolver(problem, options, writer) {

    initialize_tsppd_options();
    initialize_tsppd_constraints();
    if (valid)
        initialize_valid_inequalities();
}

void SarinTSPPDSolver::initialize_tsppd_options() {
    if (options["valid"] == "" || options["valid"] == "off")
        valid = false;
    else if (options["valid"] == "on")
        valid = true;
    else
        throw TSPPDException("valid can be either on or off");
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

void SarinTSPPDSolver::initialize_valid_inequalities() {
    // Additional valid inequalities
    for (auto pi : problem.pickup_indices()) {
        auto di = problem.successor_index(pi);
        for (auto pj : problem.pickup_indices()) {
            if (pi == pj)
                continue;

            auto dj = problem.successor_index(pj);

            model.addConstr(x[pj][di] + x[di][pj] <= y[pi][pj]);
            model.addConstr(x[pj][di] + y[di][pj] <= 1 - y[pj][pi]);

            // (+i < +j) -> (+i < -j)
            model.addConstr(y[pi][dj] >= y[pi][pj]);

            // (+i < -j) -> nothing

            // (-i < +j) -> (+i < +j) /\ (-i < -j) /\ (+i < -j)
            model.addConstr(y[pi][pj] >= y[di][pj]);
            model.addConstr(y[di][dj] >= y[di][pj]);
            model.addConstr(y[pi][dj] >= y[di][pj]);

            // (-i < -j) -> (+i < -j)
            model.addConstr(y[pi][dj] >= y[di][dj]);
        }
    }
}
