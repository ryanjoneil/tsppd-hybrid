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

#include <cmath>
#include <memory>

#include <tsppd/data/tsppd_search_statistics.h>
#include <tsppd/solver/tsp_solver.h>
#include <tsppd/solver/ruland/ruland_tsp_solver.h>
#include <tsppd/solver/ruland/callback/ruland_subtour_elimination_callback.h>
#include <tsppd/solver/ruland/callback/ruland_tsp_callback_handler.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Solver;
using namespace std;

RulandTSPSolver::RulandTSPSolver(
    const TSPPDProblem& problem,
    const map<string, string> options,
    TSPSolutionWriter& writer) :
    TSPSolver(problem, options, writer), env(), model(env), arcs(), subtour_finder(problem), callbacks() {

    // Silence output.
    model.getEnv().set(GRB_IntParam_OutputFlag, 0);

    // Lazy constraints are used for subtour elimination.
    model.set(GRB_IntParam_LazyConstraints, true);

    initialize_tsp_options();
    initialize_variables();
    initialize_two_matching_relaxation();
}

TSPPDSolution RulandTSPSolver::solve() {
    auto solver = static_cast<RulandTSPSolver*>(this);
    solver->initialize_callbacks();

    // Set time limit.
    if (time_limit > 0)
        model.set(GRB_DoubleParam_TimeLimit, time_limit / 1000.0);

    // Set solution limit.
    if (solution_limit > 0)
        model.set(GRB_IntParam_SolutionLimit, solution_limit);

    // Set thread count.
    model.set(GRB_IntParam_Threads, threads);

    RulandTSPCallbackHandler callback(solver, problem, arcs, callbacks, writer);
    model.setCallback(&callback);

    model.optimize();

    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        // Convert arc variables into boolean values for subtour finder.
        map<pair<unsigned int, unsigned int>, bool> arc_values{};
        for (auto pair : arcs)
            arc_values[pair.first] = *model.get(GRB_DoubleAttr_X, &(pair.second), 1) > 0.5;

        auto subtours = subtour_finder.subtours(arc_values);

        TSPPDSolution solution(problem, subtours[0]);
        TSPPDSearchStatistics stats(solution);
        stats.dual = ceil(model.get(GRB_DoubleAttr_ObjBound));
        writer.write(stats, true);

        return solution;
    }

    return {problem, problem.nodes};
}

void RulandTSPSolver::initialize_tsp_options() {
    if (options["sec"] == "")
        options["sec"] = "subtour";
}

void RulandTSPSolver::initialize_variables() {
    auto start_index = problem.index("+0");
    auto end_index = problem.index("-0");

    for (unsigned int from = 0; from < problem.nodes.size(); ++from) {
        for (unsigned int to = from + 1; to < problem.nodes.size(); ++to) {
            auto lb = (from == start_index && to == end_index) ? 1 : 0;
            auto var = model.addVar(lb, 1, problem.cost(from, to), GRB_BINARY);

            // Arc costs are symmetric. Recording them in both directions as the same
            // variable keeps us from having to do complicated transformations later.
            arcs[{from, to}] = var;
            arcs[{to, from}] = var;
        }
    }
}

void RulandTSPSolver::initialize_two_matching_relaxation() {
    for (unsigned int node1 = 0; node1 < problem.nodes.size(); ++node1) {
        GRBLinExpr expr = 0;
        for (unsigned int node2 = 0; node2 < problem.nodes.size(); ++node2) {
            if (node1 == node2)
                continue;

            expr += arcs[{node1, node2}];
        }

        model.addConstr(expr == 2);
    }
}

void RulandTSPSolver::initialize_callbacks() {
    callbacks.push_back(
        static_cast<shared_ptr<RulandTSPCallback>>(
            make_shared<RulandSubtourEliminationCallback>(options["sec"], problem, arcs)
        )
    );
}
