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

#ifndef TSPPD_SOLVER_RULAND_TSPPD_SOLVER_H
#define TSPPD_SOLVER_RULAND_TSPPD_SOLVER_H

#include <map>
#include <utility>

#include <gurobi_c++.h>

#include <tsppd/solver/ruland/ruland_tsp_solver.h>

namespace TSPPD {
    namespace Solver {
        // MIP TSPPD Solver based on:
        //
        // K. S. Ruland, and E. Y. Rodin. 
        // "The pickup and delivery problem: Faces and branch-and-cut algorithm." 
        // Computers & Mathematics with Applications 33, no. 12 (1997): 1-13.
        //
        // Solver Options:
        //     omc:  order matching contraints {on|off} (default=off)
        //     sec:  substour elimination constraint type
        //         - cutset:    x(delta(S)) >= 2  (default)
        //         - subtour:   sum { i,j in S } in x_{i,j} <= |S| - 1
        //         - hybrid:    subtour if |S| <= (N + 1) / 3, else cutset
        //     omc-lazy: omc lazy constraints during search {on|off} (default=off)
        class RulandTSPPDSolver : public RulandTSPSolver {
        public:
            RulandTSPPDSolver(
                const TSPPD::Data::TSPPDProblem& problem,
                const std::map<std::string, std::string> options,
                TSPPD::IO::TSPSolutionWriter& writer
            );

            virtual std::string name() const override { return "tsppd-ruland"; }

        protected:
            void initialize_tsppd_options();
            void initialize_tsppd_constraints();
            void initialize_omc_constraints();
            void initialize_omc3_constraints();
            virtual void initialize_callbacks() override;

            bool omc;
            bool omc_lazy;
       };
    }
}

#endif
