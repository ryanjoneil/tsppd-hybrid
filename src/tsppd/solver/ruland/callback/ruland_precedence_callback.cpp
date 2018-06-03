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

#include <set>

#include <tsppd/solver/ruland/callback/ruland_precedence_callback.h>
#include <tsppd/solver/ruland/callback/ruland_tsp_callback_handler.h>

using namespace TSPPD::Data;
using namespace TSPPD::Solver;
using namespace std;

RulandPrecedenceCallback::RulandPrecedenceCallback(
    const string sec_type_string,
    const TSPPDProblem& problem,
    std::map<std::pair<unsigned int, unsigned int>, GRBVar> arcs,
    const bool omc) :
    RulandSubtourEliminationCallback(sec_type_string, problem, arcs), omc(omc) { }


void RulandPrecedenceCallback::callback(GRBCallback* model, const vector<vector<unsigned int>>& subtours) {
    if (subtours.size() <= 1) {
        precedence_cut(model, subtours[0]);
        return;
    }

    for (auto tour : subtours) {
        if (omc)
            order_matching_cut(model, tour);
        cut_tour(model, tour);
    }
}

void RulandPrecedenceCallback::order_matching_cut(GRBCallback* model, const std::vector<unsigned int>& tour) {
    RulandTSPCallbackHandler* handler = static_cast<RulandTSPCallbackHandler*>(model);

    vector<int> index_pickups(problem.nodes.size(), -1);
    vector<int> index_deliveries(problem.nodes.size(), -1);
    vector<int> seen_pickups;

    for (unsigned int i = 0; i < tour.size(); ++i) {
        auto node = tour[i];

        bool is_pickup = problem.has_successor(node);
        bool is_delivery = problem.has_predecessor(node);

        if (is_pickup) {
            seen_pickups.push_back(node);
            index_pickups[node] = i;
        } else if (is_delivery)
            index_deliveries[node] = i;
    }

    for (unsigned int i = 0; i < seen_pickups.size(); ++i) {
        for (unsigned int j = i + 1 ; j < seen_pickups.size(); ++j) {
            auto p1 = seen_pickups[i];
            auto p2 = seen_pickups[j];

            auto p1_idx = index_pickups[p1];
            auto p2_idx = index_pickups[p2];

            auto d1_idx = index_deliveries[problem.successor_index(p1)];
            auto d2_idx = index_deliveries[problem.successor_index(p2)];

            if (d1_idx < 0 || d2_idx < 0 || p1_idx == 0 || p2_idx == 0)
                continue;

            if ((p1_idx - d1_idx < 0) != (p2_idx - d2_idx < 0)) {
                int start_idx = min(min(p1_idx, d1_idx), min(p2_idx, d2_idx));
                int end_idx = max(max(p1_idx, d1_idx), max(p2_idx, d2_idx));

                GRBLinExpr expr = 0;
                for (int k = start_idx; k < end_idx; ++k)
                    expr += arcs[{tour[k], tour[k+1]}];

                handler->add_lazy(expr <= end_idx - start_idx - 1);
            }
        }
    }
}

void RulandPrecedenceCallback::precedence_cut(GRBCallback* model, const vector<unsigned int>& tour) {
    RulandTSPCallbackHandler* handler = static_cast<RulandTSPCallbackHandler*>(model);

    // +0 is always first based on how we pull out the route.
    vector<unsigned int> seen_pickups(problem.nodes.size(), false);

    for (unsigned int i = 0; i < tour.size(); ++i) {
        auto node = tour[i];

        if (problem.has_successor(node))
            seen_pickups[node] = true;

        else if (problem.has_predecessor(node) && !seen_pickups[problem.predecessor_index(node)]) {
            // Precedence cuts.
            GRBLinExpr expr = 0;
            for (unsigned int node1 = 0; node1 <= i; ++node1)
                for (unsigned int node2 = i + 1; node2 < tour.size(); ++node2)
                    expr += arcs[{tour[node1], tour[node2]}];

            handler->add_lazy(expr >= 4);
        }
    }
}
