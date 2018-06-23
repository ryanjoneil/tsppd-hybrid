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

#ifndef TSPPD_SOLVER_AP_ATSPPD_CALLBACK_HANDLER_H
#define TSPPD_SOLVER_AP_ATSPPD_CALLBACK_HANDLER_H

#include <map>
#include <utility>
#include <vector>

#include <gurobi_c++.h>

#include <tsppd/data/tsppd_problem.h>
#include <tsppd/io/tsp_solution_writer.h>
#include <tsppd/solver/ap/ap_atsp_callback.h>

namespace TSPPD {
    namespace Solver {
        class APATSPPDCallback : public APATSPCallback {
        public:
            APATSPPDCallback(
                const TSPPD::Data::TSPPDProblem& problem,
                std::vector<std::vector<GRBVar>> x,
                const ATSPSECType sec,
                TSPPD::IO::TSPSolutionWriter& writer
            );

        protected:
            virtual void callback();

            std::vector<std::pair<unsigned int, unsigned int>> violations(std::vector<unsigned int> tour);
            void cut_violation(std::vector<unsigned int> tour, std::pair<unsigned int, unsigned int> index);
       };
    }
}

#endif
