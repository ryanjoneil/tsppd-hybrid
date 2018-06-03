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

#include <tsppd/solver/ruland/ruland_subtour_finder.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Solver;
using namespace std;

RulandSubtourFinder::RulandSubtourFinder(const TSPPDProblem& problem) : problem(problem) { }

vector<vector<unsigned int>> RulandSubtourFinder::subtours(
    const map<pair<unsigned int, unsigned int>, bool>& arcs) {

    vector<vector<unsigned int>> subtours{};
    auto connects = connections(arcs);

    set<unsigned int> unseen{};
    for (unsigned int node = 0; node < problem.nodes.size(); ++node)
        unseen.insert(node);

    // Start at an arbitrary node.
    auto current = *unseen.begin();
    unseen.erase(current);

    vector<unsigned int> tour = {current};

    while (!unseen.empty()) {
        // Find the next node connected to current that hasn't been seen.
        auto subtour_done = true;
        for (auto next : connects[current]) {
            if (unseen.find(next) != unseen.end()) {
                unseen.erase(next);
                subtour_done = false;
                current = next;
                tour.push_back(next);
                break;
            }
        }

        if (subtour_done) {
            subtours.push_back(tour);
            current = *unseen.begin();
            unseen.erase(current);
            tour = {current};
        }
    }

    if (tour.size() > 0)
        subtours.push_back(tour);

    // If there is only one tour, fix its direction.
    if (subtours.size() == 1) {
        vector<unsigned int> t{0};
        if (subtours[0][0] == 0 && subtours[0][1] == problem.successor_index(0))
            for (unsigned int i = subtours[0].size() - 1; i >= 1; --i)
                t.push_back(subtours[0][i]);
        else
            for (unsigned int i = 1; i < subtours[0].size(); ++i)
                t.push_back(subtours[0][i]);

        subtours[0] = t;
    }

    return subtours;
}

map<unsigned int, set<unsigned int>> RulandSubtourFinder::connections(
    const map<pair<unsigned int, unsigned int>, bool>& arcs) {

    // Create an empty set for each node's connections.
    map<unsigned int, set<unsigned int>> connects{};
    for (unsigned int node = 0; node < problem.nodes.size(); ++node)
        connects[node] = {};

    // Record every undirected connection between node pairs.
    for (unsigned int i = 0; i < problem.nodes.size(); ++i) {
        for (unsigned int j = i + 1; j < problem.nodes.size(); ++j) {
            auto arc = arcs.find({i, j});
            if (arc != arcs.end() && arc->second) {
                connects[i].insert(j);
                connects[j].insert(i);
            }
        }
    }

    return connects;
}
