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

#ifndef TSPPD_SOLVER_RULAND_PRECEDENCE_CALLBACK_H
#define TSPPD_SOLVER_RULAND_PRECEDENCE_CALLBACK_H

#include <map>
#include <vector>
#include <utility>

#include <gurobi_c++.h>

#include <tsppd/data/tsppd_problem.h>
#include <tsppd/solver/ruland/callback/ruland_subtour_elimination_callback.h>

namespace TSPPD {
    namespace Solver {
        class RulandPrecedenceCallback : public RulandSubtourEliminationCallback {
        public:
            RulandPrecedenceCallback(
                const std::string sec_type_string,
                const TSPPD::Data::TSPPDProblem& problem,
                std::map<std::pair<unsigned int, unsigned int>, GRBVar> arcs,
                const bool omc
            );

            void callback(
                GRBCallback* model,
                const std::vector<std::vector<unsigned int>>& subtours
            );

        protected:
            void order_matching_cut(GRBCallback* model, const std::vector<unsigned int>& tour);
            void precedence_cut(GRBCallback* model, const std::vector<unsigned int>& tour);

        private:
            const bool omc;
       };
    }
}

#endif
