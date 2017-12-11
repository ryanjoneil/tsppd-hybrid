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

#ifndef TSPPD_DATA_TSPPD_TREE_H
#define TSPPD_DATA_TSPPD_TREE_H

#include <memory>

#include <tsppd/data/tsppd_arc.h>
#include <tsppd/data/tsppd_problem.h>

namespace TSPPD {
    namespace Data {
        class TSPPDTree {
        public:
            TSPPDTree();

            bool contains(const unsigned int node) const;
            size_t size() const;
            void insert(const TSPPDArc& arc);
            void merge(const std::shared_ptr<TSPPDTree> other);

            static TSPPDTree minimum_spanning_tree(const TSPPDProblem& problem);

            std::set<unsigned int> nodes;
            std::vector<TSPPDArc> arcs;
            unsigned int cost;
        };

    }
}

#endif
