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

#include <tsppd/solver/sarin/sarin_atsppd_callback.h>
#include <tsppd/solver/sarin/sarin_atsppd_solver.h>
#include <tsppd/util/exception.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Solver;
using namespace TSPPD::Util;
using namespace std;

SarinATSPPDSolver::SarinATSPPDSolver(
    const TSPPDProblem& problem,
    const map<string, string> options,
    TSPSolutionWriter& writer) :
    SarinATSPSolver(problem, options, writer) {

    initialize_sarin_atsppd_options();
    initialize_sarin_atsppd_variables();

    if (valid != "none" && valid != "")
        initialize_valid_inequalities();
}

TSPPDSolution SarinATSPPDSolver::solve() {
    configure_solver();

    SarinATSPPDCallback callback(problem, x, y, sec, prec, writer);
    model.setCallback(&callback);

    model.optimize();
    return solution();
}

void SarinATSPPDSolver::initialize_sarin_atsppd_options() {
    valid = options["valid"];
    if (!(valid == "" || valid == "none" || valid == "a" || valid == "b" || valid == "all"))
        throw TSPPDException("valid can be either a, b, all, or none");

    if (options["prec"] == "" || options["prec"] == "x")
        prec = SARIN_PREC_X;
    else if (options["prec"] == "y")
        prec = SARIN_PREC_Y;
    else
        throw TSPPDException("prec can be either x or y");
}

void SarinATSPPDSolver::initialize_sarin_atsppd_variables() {
    for (auto p : problem.pickup_indices()) {
        auto d = problem.successor_index(p);

         // ~(+0 -i)
        x[start_index][d].set(GRB_DoubleAttr_UB, 0);

        // ~(+i -0)
        x[p][end_index].set(GRB_DoubleAttr_UB, 0);

        // ~(-i +i)
        x[d][p].set(GRB_DoubleAttr_UB, 0);

        // +i < -i
        y[{p,d}].set(GRB_DoubleAttr_LB, 1);
        y[{d,p}].set(GRB_DoubleAttr_UB, 0);
    }
}

void SarinATSPPDSolver::initialize_valid_inequalities() {
    // Additional valid inequalities
    for (auto pi : problem.pickup_indices()) {
        auto di = problem.successor_index(pi);
        for (auto pj : problem.pickup_indices()) {
            if (pi == pj)
                continue;

            auto dj = problem.successor_index(pj);

            if (valid == "a" || valid == "all") {
                model.addConstr(x[pj][di] + x[di][pj] <= y[{pi,pj}]);
                model.addConstr(x[pj][di] + y[{di,pj}] <= 1 - y[{pj,pi}]);
            }

            if (valid == "b" || valid == "all") {
                // (+i < +j) -> (+i < -j)
                model.addConstr(y[{pi,dj}] >= y[{pi,pj}]);

                // (+i < -j) -> nothing

                // (-i < +j) -> (+i < +j) /\ (-i < -j) /\ (+i < -j)
                model.addConstr(y[{pi,pj}] >= y[{di,pj}]);
                model.addConstr(y[{di,dj}] >= y[{di,pj}]);
                model.addConstr(y[{pi,dj}] >= y[{di,pj}]);

                // (-i < -j) -> (+i < -j)
                model.addConstr(y[{pi,dj}] >= y[{di,dj}]);
            }
        }
    }
}
