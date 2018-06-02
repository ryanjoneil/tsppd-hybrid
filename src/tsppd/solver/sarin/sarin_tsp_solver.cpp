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
#include <tsppd/solver/tsp_solver.h>
#include <tsppd/solver/sarin/sarin_tsp_solver.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Solver;
using namespace std;

SarinTSPSolver::SarinTSPSolver(
    const TSPPDProblem& problem,
    const map<string, string> options,
    TSPSolutionWriter& writer) :
    TSPSolver(problem, options, writer), env(), model(env), x(), y() {

    // Silence output.
    // model.getEnv().set(GRB_IntParam_OutputFlag, 0);

    initialize_variables();
    initialize_constraints();
}

TSPPDSolution SarinTSPSolver::solve() {
    // auto solver = static_cast<SarinTSPSolver*>(this);

    // Set time limit.
    if (time_limit > 0)
        model.set(GRB_DoubleParam_TimeLimit, time_limit / 1000.0);

    // Set solution limit.
    if (solution_limit > 0)
        model.set(GRB_IntParam_SolutionLimit, solution_limit);

    model.optimize();

    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        TSPPDSolution solution(problem, get_path());
        TSPPDSearchStatistics stats(solution);
        stats.dual = model.get(GRB_DoubleAttr_ObjBound);
        writer.write(stats, true);
        return solution;
    }

    return {problem, problem.nodes};
}

void SarinTSPSolver::initialize_variables() {
    auto start_index = problem.index("+0");
    auto end_index = problem.index("-0");

    for (unsigned int from = 0; from < problem.nodes.size(); ++from) {
        vector<GRBVar> x_i;
        vector<GRBVar> y_i;

        for (unsigned int to = 0; to < problem.nodes.size(); ++to) {
            auto fromstr = problem.nodes[from];
            auto tostr = problem.nodes[to];

            auto ub = (from == to || from == end_index || to == start_index) ? 0 : 1;
            x_i.push_back(model.addVar(0, ub, problem.cost(from, to), GRB_BINARY));

            auto lb = (from != to && to == end_index) ? 1 : 0;
            y_i.push_back(model.addVar(lb, ub, 0, GRB_BINARY));

            // TODO: BINARY? CONTINUOUS?
        }

        x.push_back(x_i);
        y.push_back(y_i);
    }
}

void SarinTSPSolver::initialize_constraints() {
    auto start_index = problem.index("+0");
    auto end_index = problem.index("-0");

    // x_ij <= y_ij
    for (unsigned int from = 0; from < problem.nodes.size(); ++from) {
        for (unsigned int to = 0; to < problem.nodes.size(); ++to) {
            if (from == to)
                continue;
            model.addConstr(x[from][to] <= y[from][to]);
        }
    }

    // sum {i != j} x_ij = 1 for all j
    for (unsigned int from = 0; from < problem.nodes.size(); ++from) {
        GRBLinExpr expr = 0;
        for (unsigned int to = 0; to < problem.nodes.size(); ++to) {
            if (from == to)
                continue;
            expr += x[from][to];
        }

        auto rhs = from == end_index ? 0 : 1;
        model.addConstr(expr == rhs);
    }

    // sum {j != i} x_ij = 1 for all i
    for (unsigned int to = 0; to < problem.nodes.size(); ++to) {
        GRBLinExpr expr = 0;
        for (unsigned int from = 0; from < problem.nodes.size(); ++from) {
            if (from == to)
                continue;
            expr += x[from][to];
        }

        auto rhs = to == start_index ? 0 : 1;
        model.addConstr(expr == rhs);
    }

    for (unsigned int i = 0; i < problem.nodes.size(); ++i)
        for (unsigned int j = i + 1; j < problem.nodes.size(); ++j)
            model.addConstr(y[i][j] + y[j][i] == 1);

    for (unsigned int i = 0; i < problem.nodes.size(); ++i) {
        for (unsigned int j = 0; j < problem.nodes.size(); ++j) {
            if (i == j)
                continue;

            for (unsigned int k = 0; k < problem.nodes.size(); ++k) {
                if (i == k || j == k)
                    continue;

                model.addConstr(y[i][j] + y[j][k] + y[k][i] <= 2);
            }
        }
    }

    // TODO: constraints and bounds to start with +0
    // TODO: lifted constraints
}

vector<unsigned int> SarinTSPSolver::get_path() {
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
