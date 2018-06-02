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
    for (auto p : problem.pickup_indices()) {
        auto d = problem.successor_index(p);
        y[p][d].set(GRB_DoubleAttr_LB, 1);
    }

    // +i < +j -> +i < -j
    for (auto p_i : problem.pickup_indices()) {
        // auto d_i = problem.successor_index(p_i);
        for (auto p_j : problem.pickup_indices()) {
            if (p_i == p_j)
                continue;
            auto d_j = problem.successor_index(p_j);
            model.addConstr(y[p_i][d_j] >= y[p_i][p_j]);
        }
    }

    // -i < -j -> +i < -j
    for (auto d_i : problem.delivery_indices()) {
        auto p_i = problem.predecessor_index(d_i);
        for (auto d_j : problem.delivery_indices()) {
            if (d_i == d_j)
                continue;
            model.addConstr(y[p_i][d_j] >= y[d_i][d_j]);
        }
    }

    // -i < +j -> +i < -j /\ -i < -j /\ +i < +j
    for (auto d_i : problem.delivery_indices()) {
        auto p_i = problem.predecessor_index(d_i);
        for (auto d_j : problem.delivery_indices()) {
            if (d_i == d_j)
                continue;
            auto p_j = problem.predecessor_index(d_j);
            
            model.addConstr(y[p_i][d_j] >= y[d_i][p_j]);
            model.addConstr(y[d_i][d_j] >= y[d_i][p_j]);
            model.addConstr(y[p_i][d_j] >= y[d_i][p_j]);
        }
    }
}

