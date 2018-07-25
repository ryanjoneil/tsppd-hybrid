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

#ifndef TSPPD_SOLVER_FOCACCI_TSPPD_SOLVER_H
#define TSPPD_SOLVER_FOCACCI_TSPPD_SOLVER_H

#include <memory>

#include <tsppd/solver/focacci/propagator/focacci_tsppd_precede_propagator.h>
#include <tsppd/solver/focacci/focacci_tsp_solver.h>
#include <tsppd/solver/focacci/focacci_tsp_space.h>

// CP TSPPD Solver based on:
//
// Filippo Focacci, Andrea Lodi, and Michela Milano.
// "A hybrid exact algorithm for the TSPTW."
// INFORMS Journal on Computing 14, no. 4 (2002): 403-417.
//
// Solver Options:
//     brancher: branching scheme {cn, regret, seq-cn} (default=regret)
//     precede:  precedence propagator type {set, cost, all} (default=set)
//     dual:     dual bounder {none, cn} (default=none)
//     omc:      order matching constraints (default=off)
//
//     search:   search engine {bab, dfs, lds} (default=bab)
//     dl:       discrepancy limit (lds only)
//
//     gist:     enables interactive search tool
//     threads:  number of threads to use in Gecode (default=1)
namespace TSPPD {
    namespace Solver {
        class FocacciTSPPDSolver : public FocacciTSPSolver {
        public:
            FocacciTSPPDSolver(
                const TSPPD::Data::TSPPDProblem& problem,
                const std::map<std::string, std::string> options,
                TSPPD::IO::TSPSolutionWriter& writer
            );

            virtual std::string name() const override { return "tsppd-cp"; }

        protected:
            void initialize_tsppd_options();
            virtual std::shared_ptr<FocacciTSPSpace> build_space() override;

            FocacciTSPPDPrecedePropagatorType precede_type;
            bool omc;
       };
    }
}

#endif
