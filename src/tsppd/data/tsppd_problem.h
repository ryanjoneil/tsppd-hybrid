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

#ifndef TSPPD_DATA_TSPPD_PROBLEM_H
#define TSPPD_DATA_TSPPD_PROBLEM_H

#include <map>
#include <queue>
#include <string>
#include <utility>
#include <vector>

#include <tsppd/data/tsppd_arc.h>

namespace TSPPD {
    namespace Data {
        enum TSPType { TSP };
        enum TSPEdgeWeightType { EUC_2D, EXPLICIT };
        enum TSPEdgeWeightFormat { NONE, LOWER_DIAG_ROW };

        class TSPPDProblem {
        public:
            TSPPDProblem();
            TSPPDProblem(
                const std::string name,
                const TSPType type,
                const std::string comment,
                const unsigned int dimension,
                const TSPEdgeWeightType edge_weight_type,
                const std::vector<std::string>& nodes,
                const std::vector<std::pair<double, double>>& coordinates,
                const std::vector<std::pair<std::string, std::string>>& pickup_delivery_pairs
            );
            TSPPDProblem(
                const std::string name,
                const TSPType type,
                const std::string comment,
                const unsigned int dimension,
                const TSPEdgeWeightType edge_weight_type,
                const std::vector<std::string>& nodes,
                const std::vector<std::pair<double, double>>& coordinates,
                const std::vector<std::pair<std::string, std::string>>& pickup_delivery_pairs,
                const std::vector<std::vector<int>>& edge_weights
            );

            void validate() const;
            void make_asymmetric(unsigned int seed);

            std::pair<double, double> coordinate(const std::string node) const;
            std::pair<double, double> coordinate(const unsigned int node_index) const;

            int cost(const std::string node1, const std::string node2) const;
            int cost(const unsigned int node_index_1, const unsigned int node_index_2) const;

            std::vector<std::string> pickups() const;
            std::vector<std::string> deliveries() const;

            std::vector<unsigned int> pickup_indices() const;
            std::vector<unsigned int> delivery_indices() const;

            unsigned int index(const std::string node) const;

            bool has_predecessor(const std::string node)        const { return has_predecessor(index(node));    }
            bool has_predecessor(const unsigned int node_index) const { return has_predecessor_vec[node_index]; }
            bool has_successor(const std::string node)          const { return has_successor(index(node));      }
            bool has_successor(const unsigned int node_index)   const { return has_successor_vec[node_index];   }

            std::string predecessor(const std::string node)        const { return predecessor(index(node));             }
            std::string predecessor(const unsigned int node_index) const { return nodes[predecessor_index(node_index)]; }
            std::string successor(const std::string node)          const { return successor(index(node));               }
            std::string successor(const unsigned int node_index)   const { return nodes[successor_index(node_index)];   }

            unsigned int predecessor_index(const unsigned int node_index) const { return predecessor_vec[node_index]; }
            unsigned int successor_index(const unsigned int node_index)   const { return successor_vec[node_index];   }

            unsigned int arcs_size(const std::string node) const;
            unsigned int arcs_size(const unsigned int node_index) const;

            std::vector<TSPPDArc> arcs(const std::string node) const;
            std::vector<TSPPDArc> arcs(const unsigned int node_index) const;

            TSPPDArc arc(const std::string node, unsigned int arc_index) const;
            TSPPDArc arc(const unsigned int node_index, unsigned int arc_index) const;

            std::string name;
            TSPType type;
            std::string comment;
            unsigned int dimension;
            TSPEdgeWeightType edge_weight_type;
            TSPEdgeWeightFormat edge_weight_format;
            std::vector<std::string> nodes;
            std::vector<std::vector<int>> edge_weights;

        private:
            std::map<std::string, unsigned int> build_node_to_index() const;

            void initialize_arc_vectors();
            void initialize_precedence(const std::vector<std::pair<std::string, std::string>>& pickup_delivery_pairs);

            std::vector<std::pair<double, double>> coordinates;

            std::vector<std::string> pickups_vec;
            std::vector<std::string> deliveries_vec;
            std::vector<unsigned int> pickup_indices_vec;
            std::vector<unsigned int> delivery_indices_vec;

            std::map<std::string, unsigned int> node_to_index;
            std::vector<std::vector<TSPPDArc>> arc_vectors;

            std::vector<bool> has_predecessor_vec;
            std::vector<bool> has_successor_vec;
            std::vector<unsigned int> predecessor_vec;
            std::vector<unsigned int> successor_vec;

            bool asymmetric = false;
       };
    }
}

#endif
