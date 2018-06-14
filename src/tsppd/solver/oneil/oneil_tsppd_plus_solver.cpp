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

#include <tsppd/solver/oneil/oneil_tsppd_plus_solver.h>
#include <tsppd/util/exception.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Solver;
using namespace TSPPD::Util;
using namespace std;

ONeilTSPPDPlusSolver::ONeilTSPPDPlusSolver(
    const TSPPDProblem& problem,
    const std::map<std::string, std::string> options,
    TSPSolutionWriter& writer) :
    ONeilTSPPDSolver(problem, options, writer),
    warm_start_solver(problem, options, writer) {

    initialize_tsppd_plus_options();

    warm_start_solver.time_limit = warm_time_limit;
    warm_start_solver.solution_limit = warm_solution_limit;

    time_limit -= warm_time_limit;
}

void ONeilTSPPDPlusSolver::initialize_tsppd_plus_options() {
    unsigned int wl = 0;
    auto wl_pair = options.find("warm-time");
    if (wl_pair != options.end())
        try {
            wl = stoi(wl_pair->second);
         } catch (exception &e) {
            throw TSPPDException("warm-time must be an integer");
         }

    warm_time_limit = wl;

    unsigned int sl = 0;
    auto sl_pair = options.find("warm-soln");
    if (sl_pair != options.end())
        try {
            sl = stoi(sl_pair->second);
         } catch (exception &e) {
            throw TSPPDException("warm-soln must be an integer");
         }

    warm_solution_limit = sl;
}

void ONeilTSPPDPlusSolver::warm_start(const TSPPDSolution& solution) {
    for (unsigned int i = 0; i < solution.order.size() - 1; ++i) {
        auto n1 = problem.index(solution.order[i]);
        auto n2 = problem.index(solution.order[i + 1]);
        x[n1][n2].set(GRB_DoubleAttr_Start, 1);
    }
}

TSPPDSolution ONeilTSPPDPlusSolver::solve() {
    warm_start(warm_start_solver.solve());
    return ONeilTSPPDSolver::solve();
}
