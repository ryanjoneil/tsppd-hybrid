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

#ifndef TSPPD_DATA_TSPPD_ARC_H
#define TSPPD_DATA_TSPPD_ARC_H

#include <limits>

namespace TSPPD {
    namespace Data {
        class TSPPDArc {
        public:
            TSPPDArc() : TSPPDArc("", "", -1, -1, std::numeric_limits<int>::max()) { }
            TSPPDArc(
                const std::string from,
                const std::string to,
                const unsigned int from_index,
                const unsigned int to_index,
                const int cost) :
                from(from), to(to), from_index(from_index), to_index(to_index), cost(cost) { }

            bool operator<(const TSPPDArc& other) const {
                return cost < other.cost;
            }

            std::string from;
            std::string to;
            unsigned int from_index;
            unsigned int to_index;
            int cost;
       };
    }
}

#endif
