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

#include <tsppd/solver/ruland/callback/ruland_subtour_elimination_callback.h>
#include <tsppd/solver/ruland/callback/ruland_tsp_callback_handler.h>
#include <tsppd/util/exception.h>

using namespace TSPPD::Data;
using namespace TSPPD::Solver;
using namespace TSPPD::Util;
using namespace std;

RulandSubtourEliminationCallback::RulandSubtourEliminationCallback(
    const string sec_type_string,
    const TSPPDProblem& problem,
    map<pair<unsigned int, unsigned int>, GRBVar> arcs) :
    RulandTSPCallback(problem, arcs), sec_type(parse_sec_type(sec_type_string)), all_nodes() {

    for (unsigned int node = 0; node < problem.nodes.size(); ++node)
        all_nodes.insert(node);
}

SECType RulandSubtourEliminationCallback::parse_sec_type(string sec_type_string) {
    if (sec_type_string == "cutset")
        return SEC_CUTSET;
    else if (sec_type_string == "subtour")
        return SEC_SUBTOUR;
    else if (sec_type_string == "hybrid")
        return SEC_HYBRID;

    throw TSPPDException("invalid sec type '" + sec_type_string + "'");
}

void RulandSubtourEliminationCallback::callback(GRBCallback* model, const vector<vector<unsigned int>>& subtours) {
    if (subtours.size() <= 1)
        return;

    for (auto tour : subtours)
        cut_tour(model, tour);
}

void RulandSubtourEliminationCallback::cut_tour(GRBCallback* model, const vector<unsigned int>& tour) {
    if (sec_type == SEC_CUTSET)
        cut_tour_cutset(model, tour);
    else if (sec_type == SEC_SUBTOUR)
        cut_tour_subtour(model, tour);
    else if (sec_type == SEC_HYBRID)
        cut_tour_hybrid(model, tour);
}

void RulandSubtourEliminationCallback::cut_tour_cutset(GRBCallback* model, const vector<unsigned int>& tour) {
    RulandTSPCallbackHandler* handler = static_cast<RulandTSPCallbackHandler*>(model);

    // Make a set with all nodes not in this subtour.
    set<unsigned int> not_in_tour(all_nodes);
    for (auto node : tour)
        not_in_tour.erase(node);

    GRBLinExpr expr = 0;
    for (auto node_in : tour)
        for (auto node_out : not_in_tour)
            expr += arcs[{node_in, node_out}];

    handler->add_lazy(expr >= 2);
}

void RulandSubtourEliminationCallback::cut_tour_subtour(GRBCallback* model, const vector<unsigned int>& tour) {
    RulandTSPCallbackHandler* handler = static_cast<RulandTSPCallbackHandler*>(model);

    GRBLinExpr expr = 0;
    for (unsigned int i = 0; i < tour.size(); ++i)
        for (unsigned int j = i + 1; j < tour.size(); ++j)
            expr += arcs[{tour[i], tour[j]}];

    handler->add_lazy(expr <= tour.size() - 1);
}

void RulandSubtourEliminationCallback::cut_tour_hybrid(GRBCallback* model, const vector<unsigned int>& tour) {
    if (tour.size() <= ((double) problem.nodes.size() + 1.0) / 3.0)
        cut_tour_subtour(model, tour);
    else
        cut_tour_cutset(model, tour);
}
