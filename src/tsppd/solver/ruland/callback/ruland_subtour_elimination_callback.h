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

#ifndef TSPPD_SOLVER_RULAND_SUBTOUR_ELIMINATION_CALLBACK_H
#define TSPPD_SOLVER_RULAND_SUBTOUR_ELIMINATION_CALLBACK_H

#include <map>
#include <vector>
#include <utility>

#include <gurobi_c++.h>

#include <tsppd/data/tsppd_problem.h>
#include <tsppd/solver/ruland/callback/ruland_tsp_callback.h>

namespace TSPPD {
    namespace Solver {
        enum SECType { SEC_CUTSET, SEC_SUBTOUR, SEC_HYBRID };

        class RulandSubtourEliminationCallback : public RulandTSPCallback {
        public:
            RulandSubtourEliminationCallback(
                const std::string sec_type_string,
                const TSPPD::Data::TSPPDProblem& problem,
                std::map<std::pair<unsigned int, unsigned int>, GRBVar> arcs
            );
            virtual ~RulandSubtourEliminationCallback() { RulandTSPCallback::~RulandTSPCallback(); }

            void callback(
                GRBCallback* model,
                const std::vector<std::vector<unsigned int>>& subtours
            );

        protected:
            SECType parse_sec_type(std::string sec_type_string);
            void cut_tour(GRBCallback* model, const std::vector<unsigned int>& tour);
            void cut_tour_cutset(GRBCallback* model, const std::vector<unsigned int>& tour);
            void cut_tour_subtour(GRBCallback* model, const std::vector<unsigned int>& tour);
            void cut_tour_hybrid(GRBCallback* model, const std::vector<unsigned int>& tour);

            const SECType sec_type;
            std::set<unsigned int> all_nodes;
       };
    }
}

#endif
