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

#include <memory>
#include <utility>
#include <vector>

#include <tsppd/data/tsppd_tree.h>

using namespace TSPPD::Data;
using namespace std;

TSPPDTree::TSPPDTree() : nodes({}), cost(0) { }

bool TSPPDTree::contains(const unsigned int node) const {
    return nodes.find(node) != nodes.end();
}

size_t TSPPDTree::size() const {
    return nodes.size();
}

void TSPPDTree::insert(const TSPPDArc& arc) {
    if (!contains(arc.from_index))
        nodes.insert(arc.from_index);

    if (!contains(arc.to_index))
        nodes.insert(arc.to_index);

    arcs.push_back(arc);
    cost += arc.cost;
}

void TSPPDTree::merge(const shared_ptr<TSPPDTree> other) {
    for (auto arc : other->arcs)
        insert(arc);
}

TSPPDTree TSPPDTree::minimum_spanning_tree(const TSPPDProblem& problem) {
    // Make an empty tree for each node and a big sorted vector of all arcs.
    vector<shared_ptr<TSPPDTree>> trees;

    vector<TSPPDArc> arcs;
    for (unsigned int i = 0; i < problem.nodes.size(); ++i) {
        auto t = make_shared<TSPPDTree>();
        t->nodes.insert(i);
        trees.push_back(t);

        for (unsigned int j = i + 1; j < problem.nodes.size(); ++j) {
            auto a = problem.arc(i, j);
            if (a.from != "")
                arcs.push_back(a);
        }
    }

    sort(arcs.begin(), arcs.end());

    TSPPDTree mst;
    for (auto arc : arcs) {
        auto tree_1 = trees[arc.from_index];
        auto tree_2 = trees[arc.to_index];
        if (tree_1 == tree_2)
            continue;

        tree_1->insert(arc);

        tree_1->merge(tree_2);
        for (auto node : tree_2->nodes)
            trees[node] = tree_1;

        if (tree_1->size() >= problem.nodes.size()) {
            mst = *tree_1;
            break;
        }
    }

    return mst;
}
