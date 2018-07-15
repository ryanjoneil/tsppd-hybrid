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

#ifndef TSPPD_SOLVER_FOCACCI_TSP_SOLVER_H
#define TSPPD_SOLVER_FOCACCI_TSP_SOLVER_H

#include <memory>

#include <tsppd/solver/tsp_solver.h>
#include <tsppd/solver/focacci/focacci_tsp_space.h>
#include <tsppd/solver/focacci/brancher/focacci_tsp_brancher.h>
#include <tsppd/solver/focacci/dual/focacci_tsp_dual.h>

// CP TSP Solver based on:
//
// Yves Caseau and Francois Laburthe.
// "Solving Small TSPs with Constraints."
// In ICLP, vol. 97, p. 104. 1997.
//
// Solver Options:
//     brancher: branching scheme {cn, regret, seq-cn} (default=regret)
//     dl:       discrepancy limit (lds only)
//     dual:     dual bounder {none, cn} (default=none)
//     gist:     enables interactive search tool (implies search=bab)
//     search:   search engine {bab, dfs, lds} (default=bab)
namespace TSPPD {
    namespace Solver {
        enum FocacciTSPSearchEngine { SEARCH_BAB, SEARCH_DFS, SEARCH_LDS };

        class FocacciTSPSolver : public TSPSolver {
        public:
            FocacciTSPSolver(
                const TSPPD::Data::TSPPDProblem& problem,
                const std::map<std::string, std::string> options,
                TSPPD::IO::TSPSolutionWriter& writer
            );

            virtual std::string name() const { return "tsp-cp"; }
            TSPPD::Data::TSPPDSolution solve();

        protected:
            void initialize_tsp_options();
            void initialize_option_brancher();
            void initialize_option_discrepancy_limit();
            void initialize_option_dual_bound();
            void initialize_option_gist();
            void initialize_option_search();

            virtual std::shared_ptr<FocacciTSPSpace> build_space();
            void check_time_limit();

            FocacciTSPBrancherType brancher_type;
            FocacciTSPDualType dual_type;
            FocacciTSPSearchEngine search_engine;

            int discrepancy_limit;
            bool gist;

            bool stopped = false;
            struct timespec start;
       };
    }
}

#endif
