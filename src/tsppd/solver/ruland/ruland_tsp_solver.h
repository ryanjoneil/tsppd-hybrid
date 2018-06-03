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

#ifndef TSPPD_SOLVER_RULAND_TSP_SOLVER_H
#define TSPPD_SOLVER_RULAND_TSP_SOLVER_H

#include <map>
#include <memory>
#include <utility>

#include <gurobi_c++.h>

#include <tsppd/solver/tsp_solver.h>
#include <tsppd/solver/ruland/ruland_subtour_finder.h>
#include <tsppd/solver/ruland/callback/ruland_tsp_callback.h>

namespace TSPPD {
    namespace Solver {
        // MIP TSP Solver based on:
        //
        // Dantzig, George, Ray Fulkerson, and Selmer Johnson. "Solution of a large-scale traveling-salesman problem." 
        // Journal of the operations research society of America 2, no. 4 (1954): 393-410.
        //
        // Solver Options:
        //     sec: subtour elimination constraint type
        //         - cutset:      x(delta(S)) >= 2
        //         - subtour:     sum { i,j in S } in x_{i,j} <= |S| - 1 (default)
        //         - hybrid:      subtour if |S| <= (N + 1) / 3, else cutset
        class RulandTSPSolver : public TSPSolver {
        public:
            RulandTSPSolver(
                const TSPPD::Data::TSPPDProblem& problem,
                const std::map<std::string, std::string> options,
                TSPPD::IO::TSPSolutionWriter& writer
            );

            std::string name() const { return "tsp-ruland"; }
            virtual TSPPD::Data::TSPPDSolution solve();

        protected:
            void initialize_tsp_options();
            void initialize_variables();
            void initialize_two_matching_relaxation();
            virtual void initialize_callbacks();

            GRBEnv env;
            GRBModel model;
            std::map<std::pair<unsigned int, unsigned int>, GRBVar> arcs;
            TSPPD::Solver::RulandSubtourFinder subtour_finder;
            std::vector<std::shared_ptr<RulandTSPCallback>> callbacks;
       };
    }
}

#endif
