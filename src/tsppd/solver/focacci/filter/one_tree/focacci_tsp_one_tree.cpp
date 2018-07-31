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
#include <iterator>
#include <limits>
#include <utility>

#include <boost/graph/kruskal_min_spanning_tree.hpp>

#include <tsppd/solver/focacci/filter/one_tree/focacci_tsp_one_tree.h>

using namespace Gecode;
using namespace TSPPD::Data;
using namespace TSPPD::Solver;
using namespace std;

OneTree::OneTree(
    ViewArray<Int::IntView>& next,
    const TSPPDProblem& problem,
    const unsigned int max_iterations) :
    OneTree(next, problem, max_iterations, nullptr) { }

OneTree::OneTree(
    ViewArray<Int::IntView>& next,
    const TSPPDProblem& problem,
    const unsigned int max_iterations,
    TSPPD::AP::PrimalDualAPSolver* ap) :
    next(next),
    problem(problem),
    max_iterations(max_iterations > 0 ? max_iterations : FOCACCI_TSP_FILTER_ONE_TREE_MAX_ITERATIONS),
    ap(ap),
    graph(next.size() - 2),
    potentials(next.size(), 0),
    edges(next.size(), set<int>()),
    start_index(problem.index("+0")),
    end_index(problem.index("-0")),
    iteration(1),
    w(0) {

    weights = boost::get(boost::edge_weight, graph);
    initialize_one_tree();
}

bool OneTree::is_done() {
    return done;
}

double OneTree::bound() {
    while (true) {
        improve();
        if (done)
            break;
    }
    return w;
}

double OneTree::improve() {
    if (iteration >= max_iterations) {
        done = true;
        return w;
    }

    // Find the max min 1-tree by updating node potentials based on violation of degree constraints.
    bool is_tour = true;

    // Compute optimal 1-tree.
    edges = vector<set<int>>(next.size(), set<int>());
    auto new_w = minimize_one_tree();

    // Update step size
    auto M = max_iterations;
    auto m = iteration;
    if (m == 1) {
        t1 = new_w / (2.0 * next.size());
        ti = t1;
    } else {
        ti = t1*(m - 1)*(2*M - 5)/(2*(M-1)) - t1*(m-2) + t1*(m-1)*(m-2)/(2*(M-1)*(M-2));
    }

    // Remove node potentials from tour.
    for (auto pi : potentials)
        new_w -= 2 * pi;

    // Update node potentials.
    for (int node = 0; node < next.size(); ++node) {
        if (node != start_index && node != end_index && edges[node].size() != 2)
            is_tour = false;
        potentials[node] += (((int) edges[node].size()) - 2) * ti;
    }

    if (is_tour || abs(new_w - w) <= EPSILON)
        done = true;
    else
        update_one_tree();

    w = new_w;
    ++iteration;
    return w;
}

bool OneTree::has_edge(int from, int to) {
    return edges[from].find(to) != edges[from].end();
}

int OneTree::marginal_cost(int from, int to) {
    vector<bool> seen(edges.size(), false);
    seen[to] = true;
    return marginal_cost(from, to, seen, to, 0);
}

void OneTree::initialize_one_tree() {
    for (int i = 0; i < next.size(); ++i) {
        // Arcs are undirected for MST.
        if (i == start_index || i == end_index)
            continue;

        for (int j = i + 1; j < next.size(); ++j) {
            if (j == start_index || j == end_index)
                continue;

            if (!(next[i].in(j) || next[j].in(i)))
                continue;

            OneTreeEdge e;
            bool inserted;
            boost::tie(e, inserted) = boost::add_edge(i, j, graph);
            weights[e] = undirected_cost(i, j);
        }
    }
}

void OneTree::update_one_tree() {
    auto es = boost::edges(graph);
    for (auto eit = es.first; eit != es.second; ++eit) {
        auto e = *eit;
        weights[e] = transformed_cost(e.m_source, e.m_target);
    }
}

double OneTree::minimize_one_tree() {
    vector<OneTreeEdge> mst;
    boost::kruskal_minimum_spanning_tree(graph, back_inserter(mst));

    // Calculate cost of the MST and pull out edges.
    double z = 0;
    for (auto e : mst) {
        edges[e.m_source].insert(e.m_target);
        edges[e.m_target].insert(e.m_source);
        z += weights[e];
    }

    // Add cheapest arc connecting +0.
    int min_p0_idx = -1;
    double min_p0 = numeric_limits<double>::max();
    for (auto to = next[start_index].min(); to <= next[start_index].max(); ++to) {
        if (!next[start_index].in(to))
            continue;

        auto c = transformed_cost(start_index, to);
        if (c < min_p0) {
            min_p0_idx = to;
            min_p0 = c;
        }
    }
    edges[start_index].insert(min_p0_idx);
    edges[min_p0_idx].insert(start_index);

    // Add cheapest arc connecting -0.
    int min_d0_idx = -1;
    double min_d0 = numeric_limits<double>::max();
    for (int from = 0; from < next.size(); ++from) {
        if (!next[from].in(end_index))
            continue;

        auto c = transformed_cost(from, end_index);
        if (c < min_d0) {
            min_d0_idx = from;
            min_d0 = c;
        }
    }
    edges[end_index].insert(min_d0_idx);
    edges[min_d0_idx].insert(end_index);

    edges[start_index].insert(end_index);
    edges[end_index].insert(start_index);

    return z + min_p0 + min_d0;
}

int OneTree::undirected_cost(int i, int j) {
    auto c_ij = next[i].in(j) ? cost(i, j) : numeric_limits<int>::max();
    auto c_ji = next[j].in(i) ? cost(j, i) : numeric_limits<int>::max();
    return min(c_ij, c_ji);
}

double OneTree::transformed_cost(int i, int j) {
    return undirected_cost(i, j) + potentials[i] + potentials[j];
}

int OneTree::marginal_cost(
    int from,
    int to,
    vector<bool> seen,
    int node,
    int max_edge_cost) {

    // Introducing a nonbasic arc into the basis would create a cycle.
    // The marginal cost of this operation is the cost of the new arc
    // minus the max cost in the cycle. This can be found using DFS.
    for (auto next : edges[node]) {
        if ((node == start_index && next == end_index) || (node == end_index && next == start_index))
            continue;

        if (seen[next])
            continue;

        // If we loop back to the from node, then compute marginal cost.
        if (next == from && node != to)
            return cost(from, to) - max_edge_cost;

        vector<bool> new_seen(seen);
        new_seen[next] = true;

        int new_max = max(max_edge_cost, undirected_cost(node, next));
        auto cost = marginal_cost(from, to, new_seen, next, new_max);
        if (cost > -1)
            return cost;
    }

    return -1;
}

int OneTree::cost(int i, int j) {
    if (ap == nullptr)
        return problem.cost(i, j);
    return ap->get_rc({i, j});
}
