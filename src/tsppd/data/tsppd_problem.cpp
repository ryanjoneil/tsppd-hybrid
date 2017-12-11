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

#include <algorithm>
#include <cmath>

#include <tsppd/data/tsppd_problem.h>
#include <tsppd/util/exception.h>

using namespace TSPPD::Data;
using namespace TSPPD::Util;
using namespace std;

TSPPDProblem::TSPPDProblem() : TSPPDProblem("", TSP, "", 0, EUC_2D, {}, {}, {}) { }

TSPPDProblem::TSPPDProblem(
    const string name,
    const TSPType type,
    const string comment,
    const unsigned int dimension,
    const TSPEdgeWeightType edge_weight_type,
    const vector<string>& nodes,
    const vector<pair<double, double>>& coordinates,
    const vector<pair<string, string>>& pickup_delivery_pairs) :
    TSPPDProblem(name, type, comment, dimension, edge_weight_type, nodes, coordinates, pickup_delivery_pairs, {}) { }

TSPPDProblem::TSPPDProblem(
    const string name,
    const TSPType type,
    const string comment,
    const unsigned int dimension,
    const TSPEdgeWeightType edge_weight_type,
    const vector<string>& nodes,
    const vector<pair<double, double>>& coordinates,
    const vector<pair<string, string>>& pickup_delivery_pairs,
    const vector<vector<int>>& edge_weights) :
    name(name),
    type(type),
    comment(comment),
    dimension(dimension),
    edge_weight_type(edge_weight_type),
    nodes(nodes),
    edge_weights(edge_weights),
    coordinates(coordinates),
    node_to_index(build_node_to_index()),
    arc_vectors(nodes.size(), vector<TSPPDArc>()),
    has_predecessor_vec(nodes.size(), false),
    has_successor_vec(nodes.size(), false),
    predecessor_vec(nodes.size(), 0),
    successor_vec(nodes.size(), 0) {

    initialize_arc_vectors();
    initialize_precedence(pickup_delivery_pairs);
}

void TSPPDProblem::validate() const {
    if (dimension < 2)
        throw TSPPDException("dimension must be >= 2");
    if (nodes.size() != dimension)
        throw TSPPDException("node count must match dimension");
    if (node_to_index.size() != dimension)
        throw TSPPDException("node indexes must be unique");
    if (coordinates.size() != dimension)
        throw TSPPDException("coordinate count must match dimension");
}

pair<double, double> TSPPDProblem::coordinate(const string node) const {
    return coordinate(index(node));
}

pair<double, double> TSPPDProblem::coordinate(const unsigned int node_index) const {
    return coordinates[node_index];
}

int TSPPDProblem::cost(const string node1, const string node2) const {
    return cost(index(node1), index(node2));
}

int TSPPDProblem::cost(const unsigned int node_index_1, const unsigned int node_index_2) const {
    if (edge_weight_type == EUC_2D) {
        auto c1 = coordinate(node_index_1);
        auto c2 = coordinate(node_index_2);
        return round(sqrt(pow(c1.first - c2.first, 2) + pow(c1.second - c2.second, 2)) + 0.5);
    }

    if (edge_weight_type == EXPLICIT) {
        return  edge_weights[max(node_index_1, node_index_2)][min(node_index_1, node_index_2)];
    }

    throw TSPPDException("invalid edge weight type");
}

unsigned int TSPPDProblem::index(const string node) const {
    auto search = node_to_index.find(node);
    if (search != node_to_index.end())
        return search->second;
    throw TSPPDException("invalid node '" + node + "'");
}

vector<string> TSPPDProblem::pickups() const {
    return pickups_vec;
}

vector<string> TSPPDProblem::deliveries() const {
    return deliveries_vec;
}

vector<unsigned int> TSPPDProblem::pickup_indices() const {
    return pickup_indices_vec;
}

vector<unsigned int> TSPPDProblem::delivery_indices() const {
    return delivery_indices_vec;
}

unsigned int TSPPDProblem::arcs_size(const std::string node) const {
    return arcs_size(index(node));
}

unsigned int TSPPDProblem::arcs_size(const unsigned int node_index) const {
    return arc_vectors[node_index].size();
}

vector<TSPPDArc> TSPPDProblem::arcs(const string node) const {
    return arcs(index(node));
}

vector<TSPPDArc> TSPPDProblem::arcs(const unsigned int node_index) const {
    return arc_vectors[node_index];
}

TSPPDArc TSPPDProblem::arc(const std::string node, unsigned int arc_index) const {
    return arc(index(node), arc_index);
}

TSPPDArc TSPPDProblem::arc(const unsigned int node_index, unsigned int arc_index) const {
    return arc_vectors[node_index][arc_index];
}

map<string, unsigned int> TSPPDProblem::build_node_to_index() const {
    map<string, unsigned int> ntoi;
    int i = 0;
    for (auto node : nodes)
        ntoi[node] = i++;
    return ntoi;
}

void TSPPDProblem::initialize_arc_vectors() {
    for (size_t i = 0; i < nodes.size(); ++i) {
        for (size_t j = i + 1; j < nodes.size(); ++j) {
            auto node1 = nodes[i];
            auto node2 = nodes[j];
            int c = cost(node1, node2);

            arc_vectors[i].push_back({node1, node2, (unsigned int) i, (unsigned int) j, c});
            arc_vectors[j].push_back({node2, node1, (unsigned int) j, (unsigned int) i, c});
        }
    }

    for (size_t i = 0; i < nodes.size(); ++i)
        sort(arc_vectors[i].begin(), arc_vectors[i].end());
}

void TSPPDProblem::initialize_precedence(const vector<pair<string, string>>& pickup_delivery_pairs) {
    for (auto p : pickup_delivery_pairs) {
        auto pickup_index = index(p.first);
        auto delivery_index = index(p.second);

        has_predecessor_vec[delivery_index] = true;
        has_successor_vec[pickup_index] = true;
        predecessor_vec[delivery_index] = pickup_index;
        successor_vec[pickup_index] = delivery_index;

        if (p.first != "+0") {
            pickups_vec.push_back(p.first);
            deliveries_vec.push_back(p.second);
            pickup_indices_vec.push_back(pickup_index);
            delivery_indices_vec.push_back(delivery_index);
        }
    }
}
