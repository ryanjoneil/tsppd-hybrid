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

#ifndef TSPPD_SOLVER_GUROBI_TSP_CALLBACK_H
#define TSPPD_SOLVER_GUROBI_TSP_CALLBACK_H

#include <map>
#include <utility>
#include <vector>

#include <gurobi_c++.h>

#include <tsppd/data/tsppd_problem.h>

namespace TSPPD {
    namespace Solver {
        class GurobiTSPCallback {
        public:
            GurobiTSPCallback(
                const TSPPD::Data::TSPPDProblem& problem,
                std::map<std::pair<unsigned int, unsigned int>, GRBVar> arcs
            ) : problem(problem), arcs(arcs) { }

            virtual void callback(
                GRBCallback* model,
                const std::vector<std::vector<unsigned int>>& subtours
            ) = 0;

        protected:
            const TSPPD::Data::TSPPDProblem& problem;
            std::map<std::pair<unsigned int, unsigned int>, GRBVar> arcs;
       };
    }
}

#endif
