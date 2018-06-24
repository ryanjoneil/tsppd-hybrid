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

#ifndef TSPPD_SOLVER_AP_ATSP_SOLVER_H
#define TSPPD_SOLVER_AP_ATSP_SOLVER_H

#include <map>

#include <gurobi_c++.h>

#include <tsppd/solver/ap/ap_atsp_callback.h>
#include <tsppd/solver/tsp_solver.h>

namespace TSPPD {
    namespace Solver {
        // MIP ATSP Solver based on an Assignment Problem relaxation
        //
        // Solver Options:
        //     sec: subtour elimination constraint type
        //         - cutset:      x(delta(S)) >= 1
        //         - subtour:     sum { i,j in S } in x_{i,j} <= |S| - 1 (default)
        class APATSPSolver : public TSPSolver {
        public:
            APATSPSolver(
                const TSPPD::Data::TSPPDProblem& problem,
                const std::map<std::string, std::string> options,
                TSPPD::IO::TSPSolutionWriter& writer
            );

            std::string name() const { return "atsp-ap"; }
            virtual TSPPD::Data::TSPPDSolution solve();

        protected:
            void initialize_options();
            void initialize_variables();
            void initialize_constraints();
            void configure_solver();

            TSPPD::Data::TSPPDSolution solution();
            std::vector<unsigned int> path();

            GRBEnv env;
            GRBModel model;
            std::vector<std::vector<GRBVar>> x;

            const unsigned int start_index;
            const unsigned int end_index;

            bool relaxed;
            ATSPSECType sec;
       };
    }
}

#endif
