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

#include <iterator>
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

    // Compute optimal 1-tree.
    vector<int> pred(next.size(), -1);
    vector<set<int>> edges(next.size(), set<int>());
    int z = one_tree(pred, edges);

    // for (unsigned int i = 0; i < edges.size(); ++i) {
    //     cout << problem.nodes[i] << ": { ";
    //     for (auto j : edges[i])
    //         cout << problem.nodes[j] << " ";
    //     cout << " }\n";
    // }

    // TODO: Perturbation of 1-tree to get HK bound

    // Objective filtering.
    GECODE_ME_CHECK(primal.gq(home, z));

    // Marginal-cost filtering,
    for (int from = 0; from < (int) next.size(); ++from) {
        for (auto to = next[from].min(); to <= next[from].max(); ++to) {
            // This only applies to nonbasic feasible arcs in the MST.
            if (!next[from].in(to) || pred[to] == from)
                continue;

            auto rc = marginal_cost(from, to, pred, edges);
            if (z + rc > primal.max())
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

int FocacciTSPPDHeldKarpPropagator::one_tree(vector<int>& pred, vector<set<int>>& edges) {
    // Construct an MST on everything but {+0,-0}.
    Graph graph(next.size() - 2);
    boost::property_map<Graph, boost::edge_weight_t>::type weights = boost::get(boost::edge_weight, graph);

    for (int i = 0; i < next.size(); ++i) {
        if (i == start_index || i == end_index)
            continue;

        for (auto j = next[i].min(); j <= next[i].max(); ++j) {
            if (!next[i].in(j) || j == (int) start_index || j == (int) end_index)
                continue;

            Edge e;
            bool inserted;
            boost::tie(e, inserted) = boost::add_edge(i, j, graph);
            weights[e] = problem.cost(i, j);
        }
    }

    vector<Edge> mst;
    boost::kruskal_minimum_spanning_tree(graph, back_inserter(mst));

    // Calculate cost of the MST and pull out the predecessor vector.
    int z = 0;
    for (auto e : mst) {
        pred[e.m_target] = e.m_source;
        edges[e.m_source].insert(e.m_target);
        edges[e.m_target].insert(e.m_source);
        z += weights[e];
    }

    // Add cheapest arc connecting +0.
    int min_p0_idx = -1;
    int min_p0 = numeric_limits<int>::max();
    for (auto to = next[start_index].min(); to <= next[start_index].max(); ++to) {
        if (!next[start_index].in(to))
            continue;

        auto c = problem.cost(start_index, to);
        if (c < min_p0) {
            min_p0_idx = to;
            min_p0 = c;
        }
    }
    pred[min_p0_idx] = start_index;
    edges[start_index].insert(min_p0_idx);
    edges[min_p0_idx].insert(start_index);

    // Add cheapest arc connecting -0.
    int min_d0_idx = -1;
    int min_d0 = numeric_limits<int>::max();
    for (int from = 0; from < next.size(); ++from) {
        if (!next[from].in(end_index))
            continue;

        auto c = problem.cost(from, end_index);
        if (c < min_d0) {
            min_d0_idx = from;
            min_d0 = c;
        }
    }
    pred[end_index] = min_d0_idx;
    edges[end_index].insert(min_d0_idx);
    edges[min_d0_idx].insert(end_index);

    return z + min_p0 + min_d0;
}

int FocacciTSPPDHeldKarpPropagator::marginal_cost(
    int from,
    int to,
    const vector<int>& pred,
    const vector<set<int>>& edges) {

    vector<bool> seen(edges.size(), false);
    seen[from] = true;
    seen[to] = true;
    return marginal_cost(from, to, pred, edges, seen, to);
}

int FocacciTSPPDHeldKarpPropagator::marginal_cost(
    int from,
    int to,
    const vector<int>& pred,
    const vector<set<int>>& edges,
    vector<bool> seen,
    int node) {

    // Introducing a nonbasic arc into the basis would create a cycle.
    // The marginal cost of this operation is the cost of the new arc
    // minus the max cost in the cycle. This can be found using DFS.
    for (auto next : edges[node]) {
        // If we loop back to the from node, then compute marginal cost.
        if (next == from) {
            if (pred[from] == node)
                return problem.cost(from, to) - problem.cost(node, from);
            return problem.cost(from, to) - problem.cost(from, node);
        }

        if (seen[next])
            continue;

        vector<bool> new_seen(seen);
        new_seen[next] = true;

        auto cost = marginal_cost(from, to, pred, edges, new_seen, next);
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
