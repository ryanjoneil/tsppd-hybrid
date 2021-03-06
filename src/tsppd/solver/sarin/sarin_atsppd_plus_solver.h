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

#ifndef TSPPD_SOLVER_SARIN_TSPPD_PLUS_SOLVER_H
#define TSPPD_SOLVER_SARIN_TSPPD_PLUS_SOLVER_H

#include <map>
#include <utility>

#include <tsppd/solver/focacci/focacci_tsppd_solver.h>
#include <tsppd/solver/sarin/sarin_atsppd_solver.h>

namespace TSPPD {
    namespace Solver {
        // MIP+CP ATSPPD Solver: Sarin MIP ATSPPD solver with time-boxed CP warm-start.
        //
        // Solver Options:
        //     relax:      relax model and add SEC and precedence as violated {on|off} (default=off)
        //     prec:       relaxed precedence form that uses either x or y variables {x|y} (default=x)
        //     sec:        relaxed SEC form that uses either x or y variables {subtour|cutset|y} (default=subtour)
        //     valid:      additional valid inequalities {a|b|all|none} (default=none)
        //     warm-time:  time limit for warm start, in milliseconds
        //     warm-soln:  solution limit for warm start
        class SarinATSPPDPlusSolver : public SarinATSPPDSolver {
        public:
            SarinATSPPDPlusSolver(
                const TSPPD::Data::TSPPDProblem& problem,
                const std::map<std::string, std::string> options,
                TSPPD::IO::TSPSolutionWriter& writer
            );

            virtual std::string name() const override { return "atsppd-sarin+"; }
            virtual TSPPD::Data::TSPPDSolution solve() override;

        protected:
            void initialize_tsppd_plus_options();
            void warm_start(const TSPPD::Data::TSPPDSolution& solution);

            FocacciTSPPDSolver warm_start_solver;

            unsigned int warm_time_limit;
            unsigned int warm_solution_limit;
        };
    }
}

#endif
