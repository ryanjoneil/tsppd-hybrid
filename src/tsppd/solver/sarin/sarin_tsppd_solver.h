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

#ifndef TSPPD_SOLVER_SARIN_TSPPD_SOLVER_H
#define TSPPD_SOLVER_SARIN_TSPPD_SOLVER_H

#include <map>

#include <tsppd/solver/sarin/sarin_tsp_solver.h>

namespace TSPPD {
    namespace Solver {
        // MIP TSPPD Solver based on:
        //
        // Subhash C. Sarin, Hanif D. Sherali, and Ajay Bhootra.
        // "New tighter polynomial length formulations for the asymmetric traveling salesman problem with
        // and without precedence constraints."
        // Operations Research Letters 33, no. 1 (2005): 62-70.
        //
        // Solver Options:
        //     relax:  relax model and add SEC and precedence as violated {on|off} (default=off)
        //     sec:    relaxed SEC form that uses either x or y variables {x|y} (default=y)
        //     valid:  additional valid inequalities {a|b|all|none} (default=none)
        class SarinTSPPDSolver : public SarinTSPSolver {
        public:
            SarinTSPPDSolver(
                const TSPPD::Data::TSPPDProblem& problem,
                const std::map<std::string, std::string> options,
                TSPPD::IO::TSPSolutionWriter& writer
            );

            std::string name() const { return "tsppd-sarin"; }

        protected:
            void initialize_tsppd_options();
            void initialize_tsppd_constraints();
            void initialize_valid_inequalities();

            std::string valid;
       };
    }
}

#endif
