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

#include <set>

#include <tsppd/data/tsppd_solution.h>

using namespace TSPPD::Data;
using namespace std;

TSPPDSolution::TSPPDSolution(const TSPPDProblem& problem, const vector<string> order) :
    problem(problem), order(order), cost(compute_cost()) { }

TSPPDSolution::TSPPDSolution(const TSPPD::Data::TSPPDProblem& problem, const std::vector<unsigned int> index_order) :
    problem(problem), order(string_order(index_order)), cost(compute_cost()) { }

vector<string> TSPPDSolution::string_order(vector<unsigned int> index_order) {
    vector<string> tour;
    for (auto i : index_order)
        tour.push_back(problem.nodes[i]);
    return tour;
}

bool TSPPDSolution::feasible() const {
    // Check size and check precedence.
    if (order.size() != problem.nodes.size())
        return false;

    set<string> seen;
    for (auto node : order) {
        seen.insert(node);
        if (problem.has_predecessor(node) && seen.find(problem.predecessor(node)) == seen.end())
            return false;
    }

    return true;
}

int TSPPDSolution::compute_cost() {
    // Compute total cost of the tour.
    int tour_cost = 0;
    for (size_t i = 0; i < order.size() - 1; ++i)
        tour_cost += problem.cost(order[i], order[i+1]);

    // Make sure to include the final leg back to the origin.
    tour_cost += problem.cost(order[order.size() - 1], order[0]);
    return tour_cost;
}
