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

#include <cstdio>
#include <vector>

#include <tsppd/data/tsppd_search_statistics.h>
#include <tsppd/solver/oneil/oneil_tsppd_callback.h>
#include <tsppd/solver/oneil/oneil_tsppd_solver.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Solver;
using namespace std;

ONeilTSPPDSolver::ONeilTSPPDSolver(
    const TSPPDProblem& problem,
    const map<string, string> options,
    TSPSolutionWriter& writer) :
    TSPSolver(problem, options, writer),
    env(),
    model(env),
    x(),
    y(),
    start_index(problem.index("+0")),
    end_index(problem.index("-0")) {

    // Silence output.
    model.getEnv().set(GRB_IntParam_OutputFlag, 0);

    initialize_variables();
    initialize_assignment_problem_constraints();
    initialize_subtour_and_precedence_constraints();
}

TSPPDSolution ONeilTSPPDSolver::solve() {
    // Set time limit.
    if (time_limit > 0)
        model.set(GRB_DoubleParam_TimeLimit, time_limit / 1000.0);

    // Set solution limit.
    if (solution_limit > 0)
        model.set(GRB_IntParam_SolutionLimit, solution_limit);

    // Set thread count.
    model.set(GRB_IntParam_Threads, threads);

    ONeilTSPPDCallback callback(problem, x, writer);
    model.setCallback(&callback);

    model.optimize();

    // cout << endl;

    // cout << "y = " << endl;
    // for (unsigned int i = 0; i < problem.nodes.size(); ++i) {
    //     cout << "[" << problem.nodes[i] << "] ";
    //     for (unsigned int j = 0; j < problem.nodes.size(); ++j)
    //         if (*model.get(GRB_DoubleAttr_X, &(y[i][j]), 1) > 0.5)
    //             cout << problem.nodes[j] << " ";
    //     cout << endl;
    // }

    // cout << endl;

    // cout << "x = " << endl;
    // for (unsigned int i = 0; i < problem.nodes.size(); ++i) {
    //     cout << "[" << problem.nodes[i] << "] ";
    //     for (unsigned int j = 0; j < problem.nodes.size(); ++j)
    //         if (*model.get(GRB_DoubleAttr_X, &(x[i][j]), 1) > 0.5)
    //             cout << problem.nodes[j] << " ";
    //     cout << endl;
    // }

    // cout << endl;


    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        TSPPDSolution solution(problem, get_path());
        TSPPDSearchStatistics stats(solution);
        stats.dual = model.get(GRB_DoubleAttr_ObjBound);
        writer.write(stats, true);
        return solution;
    }

    return {problem, problem.nodes};
}

void ONeilTSPPDSolver::initialize_variables() {
    for (unsigned int from = 0; from < problem.nodes.size(); ++from) {
        vector<GRBVar> x_i;
        vector<GRBVar> y_i;

        for (unsigned int to = 0; to < problem.nodes.size(); ++to) {
            auto fromstr = problem.nodes[from];
            auto tostr = problem.nodes[to];

            auto ub_y = (from == to || from == end_index || to == start_index) ? 0 : 1;
            auto lb_y = (from != to && (from == start_index || to == end_index)) ? 1 : 0;
            y_i.push_back(model.addVar(lb_y, ub_y, 0, GRB_BINARY,
                string("y[" + problem.nodes[from] + " " + problem.nodes[to] + "]")));

            // +0 can only connect to a pickup, only deliveries can connect to -0
            auto ub_x = ub_y;
            if (from == start_index && !problem.has_successor(to))
                ub_x = 0;
            else if (to == end_index && !problem.has_predecessor(from))
                ub_x = 0;

            x_i.push_back(model.addVar(0, ub_x, problem.cost(from, to), GRB_BINARY,
                string("x[" + problem.nodes[from] + " " + problem.nodes[to] + "]")));

        }

        x.push_back(x_i);
        y.push_back(y_i);
    }
}

void ONeilTSPPDSolver::initialize_assignment_problem_constraints() {
    // x_ij <= y_ij for all i,j = 2,...,n, i != j
    for (unsigned int from = 0; from < problem.nodes.size(); ++from) {
        if (from == start_index)
            continue;

        for (unsigned int to = 0; to < problem.nodes.size(); ++to) {
            if (to == from || to == start_index)
                continue;

            model.addConstr(x[from][to] <= y[from][to]);
        }
    }

    // sum {j != i} x_ij = 1 for all i = 2,...,n
    for (unsigned int to = 0; to < problem.nodes.size(); ++to) {
        if (to == start_index)
            continue;

        GRBLinExpr expr = 0;
        for (unsigned int from = 0; from < problem.nodes.size(); ++from) {
            if (from == to)
                continue;
            expr += x[from][to];
        }

        // +0 is preceded by nothing
        auto rhs = to == start_index ? 0 : 1;
        model.addConstr(expr == rhs);
    }

    // sum {i != j} x_ij = 1 for all j = 2,...,n
    for (unsigned int from = 0; from < problem.nodes.size(); ++from) {
        if (from == start_index)
            continue;

        GRBLinExpr expr = 0;
        for (unsigned int to = 0; to < problem.nodes.size(); ++to) {
            if (from == to)
                continue;

            expr += x[from][to];
        }

        // -0 precedes nothing
        auto rhs = from == end_index ? 0 : 1;
        model.addConstr(expr == rhs);
    }
}

void ONeilTSPPDSolver::initialize_subtour_and_precedence_constraints() {
    // +i < -i
    for (auto p : problem.pickup_indices()) {
        auto d = problem.successor_index(p);
        y[p][d].set(GRB_DoubleAttr_LB, 1);
    }

    // y_ij + y_ji = 1
    for (unsigned int i = 0; i < problem.nodes.size(); ++i) {
        for (unsigned int j = i + 1; j < problem.nodes.size(); ++j) {
            if (j == start_index)
                continue;

            model.addConstr(y[i][j] + y[j][i] == 1);
        }
    }

    // +i < +j -> +i < -j
    for (auto p_i : problem.pickup_indices()) {
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

    // y[i][j] /\ y[j][k] -> y[i][k]
    for (unsigned int i = 0; i < problem.nodes.size(); ++i) {
        if (i == start_index || i == end_index)
            continue;

        for (unsigned int j = 0; j < problem.nodes.size(); ++j) {
            if (i == start_index || i == end_index || j == i)
                continue;

            for (unsigned int k = 0; k < problem.nodes.size(); ++k) {
                if (i == start_index || i == end_index || k == i || k == j)
                    continue;

                model.addConstr(y[i][k] + x[j][i] >= y[i][j] + y[j][k] - 1);
            }
        }
    }    
}

vector<unsigned int> ONeilTSPPDSolver::get_path() {
    vector<unsigned int> path{0};
    unsigned int from = 0;
    for (unsigned int i = 0; i < problem.nodes.size() - 1; ++i) {
        for (unsigned int to = 0; to < problem.nodes.size(); ++to) {
            if (*model.get(GRB_DoubleAttr_X, &(x[from][to]), 1) > 0.5) {
                path.push_back(to);
                from = to;
                break;
            }
        }
    }
    return path;
}
