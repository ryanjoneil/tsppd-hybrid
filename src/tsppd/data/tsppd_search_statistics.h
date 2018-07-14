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

#ifndef TSPPD_DATA_TSPPD_SEARCH_STATISTICS_H
#define TSPPD_DATA_TSPPD_SEARCH_STATISTICS_H

#include <vector>

#include <tsppd/data/tsppd_solution.h>

namespace TSPPD {
    namespace Data {
        class TSPPDSearchStatistics {
        public:
            TSPPDSearchStatistics() { }
            TSPPDSearchStatistics(const TSPPD::Data::TSPPDSolution& solution) :
                primal(solution.cost), tour(solution.order) { }

            // Solution data
            bool has_primal() const { return primal >= 0; }
            bool has_dual() const { return dual >= 0; }
            bool has_tour() const { return tour.size() > 0; }
            bool is_optimal() const { return optimal; }

            // Search information
            bool has_nodes() const { return nodes >= 0; }
            bool has_fails() const { return fails >= 0; }
            bool has_depth() const { return depth >= 0; }

            int primal = -1;
            int dual = -1;
            std::vector<std::string> tour = {};
            bool optimal = false;

            int nodes = -1;
            int fails = -1;
            int depth = -1;

        };
    }
}

#endif
