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

#ifndef TSPPD_SOLVER_AP_ATSP_CALLBACK_HANDLER_H
#define TSPPD_SOLVER_AP_ATSP_CALLBACK_HANDLER_H

#include <utility>
#include <vector>

#include <gurobi_c++.h>

#include <tsppd/data/tsppd_problem.h>
#include <tsppd/io/tsp_solution_writer.h>

namespace TSPPD {
    namespace Solver {
        enum ATSPSECType { ATSP_SEC_CUTSET, ATSP_SEC_SUBTOUR, ATSP_SEC_OTHER };

        class APATSPCallback : public GRBCallback {
        public:
            APATSPCallback(
                const TSPPD::Data::TSPPDProblem& problem,
                std::vector<std::vector<GRBVar>> x,
                const ATSPSECType sec_type,
                TSPPD::IO::TSPSolutionWriter& writer
            );

        protected:
            virtual void callback();
            void log_mip();
            void log_mipsol(std::vector<unsigned int>& tour);

            std::vector<std::vector<unsigned int>> subtours();
            virtual void cut_subtour(const std::vector<unsigned int>& subtour);
            void cut_subtour_cutset(const std::vector<unsigned int>& subtour);
            void cut_subtour_subtour(const std::vector<unsigned int>& subtour);

            const TSPPD::Data::TSPPDProblem& problem;
            std::vector<std::vector<GRBVar>> x;
            const ATSPSECType sec_type;
            TSPPD::IO::TSPSolutionWriter writer;
       };
    }
}

#endif
