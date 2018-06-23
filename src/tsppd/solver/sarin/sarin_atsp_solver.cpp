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

#include <tsppd/solver/sarin/sarin_atsp_callback.h>
#include <tsppd/solver/sarin/sarin_atsp_solver.h>
#include <tsppd/util/exception.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Solver;
using namespace TSPPD::Util;
using namespace std;

SarinATSPSolver::SarinATSPSolver(
    const TSPPDProblem& problem,
    const map<string, string> options,
    TSPSolutionWriter& writer) :
    APATSPSolver(problem, options, writer),
    y() {

    initialize_sarin_options();
    initialize_sarin_variables();
    initialize_sarin_constraints();
}

TSPPDSolution SarinATSPSolver::solve() {
    configure_solver();

    SarinATSPCallback callback(problem, x, y, sec, writer);
    model.setCallback(&callback);

    model.optimize();
    return solution();
}

void SarinATSPSolver::initialize_sarin_options() {
    if (options["relax"] == "" || options["relax"] == "off")
        relaxed = false;
    else if (options["relax"] == "on")
        relaxed = true;
    else
        throw TSPPDException("relax can be either on or off");

    if (options["sec"] == "" || options["sec"] == "subtour")
        sec = ATSP_SEC_SUBTOUR;
    else if (options["sec"] == "cutset")
        sec = ATSP_SEC_CUTSET;
    else if (options["sec"] == "y")
        sec = ATSP_SEC_OTHER;
    else
        throw TSPPDException("sec can be either subtour, cutset, or y");
}

void SarinATSPSolver::initialize_sarin_variables() {
    for (unsigned int from = 0; from < problem.nodes.size(); ++from) {
        if (from == start_index || from == end_index)
            continue;

        // if (problem.has_predecessor(from))
        //     x[start_index][from].set(GRB_DoubleAttr_UB, 0);
        // else if (problem.has_successor(from))
        //     x[from][end_index].set(GRB_DoubleAttr_UB, 0);

        for (unsigned int to = 0; to < problem.nodes.size(); ++to) {
            if (to == start_index || to == end_index || to == from)
                continue;

            auto ub = 1;
            // if (problem.has_successor(from) && problem.successor_index(from) == to)
            //     ub = 0;

            y[{from,to}] = model.addVar(0, ub, 0, GRB_CONTINUOUS);
        }
    }
}

void SarinATSPSolver::initialize_sarin_constraints() {
    // x_ij <= y_ij for all i,j
    for (auto& yi : y) {
        auto ind = yi.first;
        auto v = yi.second;

        auto from = ind.first;
        auto to = ind.second;
        model.addConstr(x[from][to] <= v);

        if (from < to)
            model.addConstr(y[{from,to}] + y[{to,from}] == 1);
    }

    if (!relaxed) {
        // (y_ij + x_ji) + y_jk + y_ki <= 2 for all i,j,k = 2,...,n, i != j != k
        for (unsigned int i = 0; i < problem.nodes.size(); ++i) {
            if (i == start_index || i == end_index)
                continue;

            for (unsigned int j = 0; j < problem.nodes.size(); ++j) {
                if (i == j || j == start_index || j == end_index)
                    continue;

                for (unsigned int k = 0; k < problem.nodes.size(); ++k) {
                    if (i == k || j == k || k == start_index || k == end_index)
                        continue;

                    model.addConstr(y[{i,j}] + x[j][i] + y[{j,k}] + y[{k,i}] <= 2);
                }
            }
        }
    }
}
