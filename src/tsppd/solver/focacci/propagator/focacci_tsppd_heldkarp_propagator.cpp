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

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>

#include <tsppd/solver/focacci/propagator/focacci_tsppd_heldkarp_propagator.h>

using namespace Gecode;
using namespace TSPPD::AP;
using namespace TSPPD::Data;
using namespace TSPPD::Solver;
using namespace std;

typedef boost::adjacency_list<
    boost::vecS,
    boost::vecS,
    boost::directedS,
    boost::no_property,
    boost::property<boost::edge_weight_t, int>
> Graph;
typedef boost::graph_traits <Graph>::edge_descriptor Edge;

FocacciTSPPDHeldKarpPropagator::FocacciTSPPDHeldKarpPropagator(
    Home home,
    ViewArray<Int::IntView>& next,
    Int::IntView& primal,
    const TSPPDProblem& problem) :
    Propagator(home),
    next(next),
    primal(primal),
    problem(problem),
    start_index(problem.index("+0")),
    end_index(problem.index("-0")) {

    next.subscribe(home, *this, Int::PC_INT_DOM);
    home.notice(*this, AP_DISPOSE);
}

FocacciTSPPDHeldKarpPropagator::FocacciTSPPDHeldKarpPropagator(Space& home, FocacciTSPPDHeldKarpPropagator& p) :
    Propagator(home, p),
    next(p.next),
    primal(p.primal),
    problem(p.problem),
    start_index(p.start_index),
    end_index(p.end_index) {

    next.update(home, p.next);
    primal.update(home, p.primal);
}

Propagator* FocacciTSPPDHeldKarpPropagator::copy(Space& home) {
    return new (home) FocacciTSPPDHeldKarpPropagator(home, *this);
}

size_t FocacciTSPPDHeldKarpPropagator::dispose(Space& home) {
    home.ignore(*this, AP_DISPOSE);
    next.cancel(home, *this, Int::PC_INT_DOM);
    (void) Propagator::dispose(home);
    return sizeof(*this);
}

PropCost FocacciTSPPDHeldKarpPropagator::cost(const Space& home, const ModEventDelta& med) const {
    return PropCost::crazy(PropCost::HI, next.size());
}

void FocacciTSPPDHeldKarpPropagator::reschedule(Space& home) {
    next.reschedule(home, *this, Int::PC_INT_DOM);
}

ExecStatus FocacciTSPPDHeldKarpPropagator::propagate(Space& home, const ModEventDelta& med) {
    if (primal.assigned() || next.assigned())
        return home.ES_SUBSUMED(*this);

    vector<double> potentials(next.size(), 0);
    vector<set<int>> edges;
    double w;
    double last_w = 0;

    for (unsigned int iter = 0; iter < MAX_ITERATIONS; ++iter) {
        bool is_tour = true;

        // Compute optimal 1-tree.
        edges = vector<set<int>>(next.size(), set<int>());
        w = one_tree(potentials, edges);

        // Remove node potentials from tour.
        for (auto pi : potentials)
            w -= 2 * pi;

        // Update node potentials.
        for (int node = 0; node < next.size(); ++node) {
            if (node != start_index && node != end_index && edges[node].size() != 2)
                is_tour = false;
            potentials[node] += (((int) edges[node].size()) - 2) * STEP_SIZE;
        }

        // cout << "[" << iter << "] 1-tree=" << w << " primal=" << primal << " t=" << t << endl;
        if (is_tour || w <= last_w + 0.1)
            break;
        else
            last_w = w;
    }

    // cout << "1-tree=" << w << " primal=" << primal << " t=" << t << endl;
    // for (unsigned int i = 0; i < edges.size(); ++i) {
    //     cout << problem.nodes[i] << ": { ";
    //     for (auto j : edges[i])
    //         cout << problem.nodes[j] << " ";
    //     cout << " }\n";
    // }

    // Perturbation of 1-tree to get HK bound

    // Objective filtering.
    GECODE_ME_CHECK(primal.gq(home, (int) ceil(w)));

    // Marginal-cost filtering.
    // TODO: what do we do with node potentials?
    for (int from = 0; from < (int) next.size(); ++from) {
        for (auto to = next[from].min(); to <= next[from].max(); ++to) {
            // This only applies to nonbasic feasible arcs in the MST.
            if (!(next[from].in(to) && edges[from].find(to) == edges[from].end()))
                continue;

            if (w + marginal_cost(from, to, edges) > primal.max())
                GECODE_ME_CHECK(next[from].nq(home, to));
        }
    }

    return ES_FIX;
}

ExecStatus FocacciTSPPDHeldKarpPropagator::post(
    Home home,
    ViewArray<Int::IntView>& next,
    Int::IntView& primal,
    const TSPPDProblem& problem) {

    if (!primal.assigned() && !next.assigned())
        (void) new (home) FocacciTSPPDHeldKarpPropagator(home, next, primal, problem);
    return ES_OK;
}

double FocacciTSPPDHeldKarpPropagator::one_tree(const vector<double>& potentials, vector<set<int>>& edges) {
    // Construct an MST on everything but {+0,-0}.
    Graph graph(next.size() - 2);
    boost::property_map<Graph, boost::edge_weight_t>::type weights = boost::get(boost::edge_weight, graph);

    for (int i = 0; i < next.size(); ++i) {
        // Arcs are undirected for MST.
        if (i == start_index || i == end_index)
            continue;

        for (int j = i + 1; j < next.size(); ++j) {
            if (j == start_index || j == end_index)
                continue;

            if (!(next[i].in(j) || next[j].in(i)))
                continue;

            Edge e;
            bool inserted;
            boost::tie(e, inserted) = boost::add_edge(i, j, graph);
            weights[e] = transformed_cost(i, j, potentials);
        }
    }

    vector<Edge> mst;
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

        auto c = transformed_cost(start_index, to, potentials);
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

        auto c = transformed_cost(from, end_index, potentials);
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

int FocacciTSPPDHeldKarpPropagator::undirected_cost(int i, int j) {
    auto c_ij = next[i].in(j) ? problem.cost(i, j) : numeric_limits<int>::max();
    auto c_ji = next[j].in(i) ? problem.cost(j, i) : numeric_limits<int>::max();
    return min(c_ij, c_ji);
}

double FocacciTSPPDHeldKarpPropagator::transformed_cost(int i, int j, const std::vector<double> potentials) {
    return undirected_cost(i, j) + potentials[i] + potentials[j];
}

// TODO: should this operate on original or transformed costs?
int FocacciTSPPDHeldKarpPropagator::marginal_cost(
    int from,
    int to,
    const vector<set<int>>& edges) {

    vector<bool> seen(edges.size(), false);
    seen[to] = true;
    return marginal_cost(from, to, edges, seen, to, 0);
}

int FocacciTSPPDHeldKarpPropagator::marginal_cost(
    int from,
    int to,
    const vector<set<int>>& edges,
    vector<bool> seen,
    int node,
    int max_edge_cost) {

    // Introducing a nonbasic arc into the basis would create a cycle.
    // The marginal cost of this operation is the cost of the new arc
    // minus the max cost in the cycle. This can be found using DFS.
    for (auto next : edges[node]) {
        if (seen[next])
            continue;

        int new_max = max(new_max, undirected_cost(from, node));

        // If we loop back to the from node, then compute marginal cost.
        if (next == from && node != to)
            return problem.cost(from, to) - new_max;

        vector<bool> new_seen(seen);
        new_seen[next] = true;

        auto cost = marginal_cost(from, to, edges, new_seen, next, new_max);
        if (cost > -1)
            return cost;
    }

    return -1;
}

void TSPPD::Solver::tsppd_heldkarp(
    Home home,
    IntVarArray& next,
    IntVar& primal,
    const TSPPDProblem& problem) {

    GECODE_POST;

    IntVarArgs next_args(next);
    ViewArray<Int::IntView> next_view(home, next_args);

    Int::IntView primal_view(primal);

    GECODE_ES_FAIL(FocacciTSPPDHeldKarpPropagator::post(home, next_view, primal_view, problem));
}
