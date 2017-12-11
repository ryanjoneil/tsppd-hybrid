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

#include <tsppd/solver/enumerative/enumerative_tsppd_solver.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Solver;
using namespace std;

// Enumerative TSP Solver
EnumerativeTSPPDSolver::EnumerativeTSPPDSolver(
    const TSPPDProblem& problem,
    const map<string, string> options,
    TSPSolutionWriter& writer) :
    EnumerativeTSPSolver(problem, options, writer) { }

bool EnumerativeTSPPDSolver::feasible(TSPPDArc next) {
    if (problem.has_predecessor(next.to_index)) {
        auto pickup = problem.predecessor_index(next.to_index);
        if (!in_tour[pickup]) {
            return false;
        }

        // End node can only be feasible at end
        if (pickup == 0 && current_tour.size() < problem.nodes.size() - 1) {
            return false;
        }
    }

    return EnumerativeTSPSolver::feasible(next);
}
