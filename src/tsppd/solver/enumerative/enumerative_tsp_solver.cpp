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

#include <limits>

#include <tsppd/data/tsppd_search_statistics.h>
#include <tsppd/solver/enumerative/enumerative_tsp_solver.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Solver;
using namespace std;

// Enumerative TSP Solver
EnumerativeTSPSolver::EnumerativeTSPSolver(
    const TSPPDProblem& problem,
    const map<string, string> options,
    TSPSolutionWriter& writer) :
    TSPSolver(problem, options, writer),
    arcs(),
    in_tour(problem.nodes.size(), false),
    current_tour(),
    start(clock()) {

    initialize_search();
}

TSPPDSolution EnumerativeTSPSolver::solve() {
    find_best();

    TSPPDSolution solution(problem, best_tour);
    TSPPDSearchStatistics stats(solution);
    writer.write(stats, true);

    return solution;
}

void EnumerativeTSPSolver::initialize_search() {
    for (unsigned int i = 0; i < problem.nodes.size(); ++i)
        arcs.push_back(problem.arcs(i));
    current_tour.push_back(0);
    in_tour[0] = true;
    current_cost = 0;
    best_cost = numeric_limits<int>::max();
}

void EnumerativeTSPSolver::find_best() {
    if (stopped)
        return;

    auto current = current_tour[current_tour.size() - 1];

    for (unsigned int index = 0; index < arcs[current].size(); ++index) {
        // Time limits
        if (time_limit > 0 && ((clock() - start) * 1000.0 / CLOCKS_PER_SEC) >= time_limit) {
            stopped = true;
            return;
        }

        auto next = arcs[current][index];
        if (!feasible(next))
            continue;

        // Continue exploring this part of the tree.
        current_tour.push_back(next.to_index);
        current_cost += next.cost;
        in_tour[next.to_index] = true;

        if (current_tour.size() > problem.nodes.size() - 1) {
            // Add in final arc.
            int final_arc_cost = problem.cost(0, current_tour[current_tour.size() - 1]);
            current_cost += final_arc_cost;

            if (current_cost < best_cost) {
                best_tour = current_tour;
                best_cost = current_cost;

                TSPPDSolution solution(problem, current_tour);
                TSPPDSearchStatistics stats(solution);
                writer.write(solution);

                // If we have a solution limit, respect it.
                if (solution_limit > 0) {
                    if (--solution_limit == 0) {
                        stopped = true;
                        return;
                    }
                }
            }

            current_cost -= final_arc_cost;

        } else {
            find_best();
        }

        current_tour.pop_back();
        current_cost -= next.cost;
        in_tour[next.to_index] = false;
    }
}

bool EnumerativeTSPSolver::feasible(TSPPDArc next) {
    // Node is already in the tour.
    if (in_tour[next.to_index])
        return false;

    // Fathom sections of the search tree that give worse results.
    if (current_cost + next.cost >= best_cost)
        return false;

    return true;
}
