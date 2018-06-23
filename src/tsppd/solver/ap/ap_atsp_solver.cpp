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

#include <tsppd/data/tsppd_search_statistics.h>
#include <tsppd/solver/ap/ap_atsp_callback.h>
#include <tsppd/solver/ap/ap_atsp_solver.h>
#include <tsppd/util/exception.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Solver;
using namespace TSPPD::Util;
using namespace std;

APATSPSolver::APATSPSolver(
    const TSPPDProblem& problem,
    const map<string, string> options,
    TSPSolutionWriter& writer) :
    TSPSolver(problem, options, writer),
    env(),
    model(env),
    x(),
    start_index(problem.index("+0")),
    end_index(problem.index("-0")) {

    initialize_options();
    initialize_variables();
    initialize_constraints();
}

TSPPDSolution APATSPSolver::solve() {
    configure_solver();

    auto sec_type = sec == "cutset" ? AP_ATSP_SEC_CUTSET : AP_ATSP_SEC_SUBTOUR;
    APATSPCallback callback(problem, x, sec_type, writer);
    model.setCallback(&callback);

    model.optimize();
    return solution();
}

void APATSPSolver::initialize_options() {
    relaxed = true;

    if (options["sec"] == "" || options["sec"] == "subtour")
        sec = "subtour";
    else if (options["sec"] == "cutset")
        sec = "cutset";
    else
        throw TSPPDException("sec can be either subtour or cutset");
}

void APATSPSolver::initialize_variables() {
    for (unsigned int from = 0; from < problem.nodes.size(); ++from) {
        vector<GRBVar> x_i;

        for (unsigned int to = 0; to < problem.nodes.size(); ++to) {
            auto lb = (from == end_index && to == start_index) ? 1 : 0;
            auto ub = (from == to) ? 0 : 1;
            x_i.push_back(model.addVar(lb, ub, problem.cost(from, to), GRB_BINARY));
        }

        x.push_back(x_i);
    }
}

void APATSPSolver::initialize_constraints() {
    // sum {j != i} x_ij = 1 for all i
    for (unsigned int to = 0; to < problem.nodes.size(); ++to) {
        GRBLinExpr expr = 0;
        for (unsigned int from = 0; from < problem.nodes.size(); ++from) {
            if (from == to)
                continue;
            expr += x[from][to];
        }
        model.addConstr(expr == 1);
    }

    // sum {i != j} x_ij = 1 for all j
    for (unsigned int from = 0; from < problem.nodes.size(); ++from) {
        GRBLinExpr expr = 0;
        for (unsigned int to = 0; to < problem.nodes.size(); ++to) {
            if (from == to)
                continue;
            expr += x[from][to];
        }
        model.addConstr(expr == 1);
    }
}

void APATSPSolver::configure_solver() {
    // Silence output.
    model.getEnv().set(GRB_IntParam_OutputFlag, 0);

    // Lazy constraints are used for subtour elimination.
    if (relaxed)
        model.set(GRB_IntParam_LazyConstraints, true);

    // Set time limit.
    if (time_limit > 0)
        model.set(GRB_DoubleParam_TimeLimit, time_limit / 1000.0);

    // Set solution limit.
    if (solution_limit > 0)
        model.set(GRB_IntParam_SolutionLimit, solution_limit);

    // Set thread count.
    model.set(GRB_IntParam_Threads, threads);
}

TSPPDSolution APATSPSolver::solution() {
    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        TSPPDSolution solution(problem, path());
        TSPPDSearchStatistics stats(solution);
        stats.dual = model.get(GRB_DoubleAttr_ObjBound);
        writer.write(stats, true);
        return solution;
    }

    return {problem, problem.nodes};
}

vector<unsigned int> APATSPSolver::path() {
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
