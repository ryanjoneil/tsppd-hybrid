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

#include <vector>

#include <tsppd/data/tsppd_search_statistics.h>
#include <tsppd/data/tsppd_solution.h>
#include <tsppd/solver/ap/ap_atsppd_callback.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Solver;
using namespace std;

APATSPPDCallback::APATSPPDCallback(
    const TSPPDProblem& problem,
    vector<vector<GRBVar>> x,
    const ATSPSECType sec,
    TSPSolutionWriter& writer) :
    APATSPCallback(problem, x, sec, writer) { }

void APATSPPDCallback::callback() {
    if (where == GRB_CB_MIP)
        log_mip();

    else if (where == GRB_CB_MIPSOL) {
        auto s = subtours();

        if (s.size() > 1) {
            for (auto subtour : s)
                cut_subtour(subtour);
        } else {
            auto v = violations(s.front());
            if (v.size() > 0)
                for (auto vi : v)
                    cut_violation(s.front(), vi);
            else
                log_mipsol(s.front());
        }
    }
}

vector<pair<unsigned int, unsigned int>>APATSPPDCallback::violations(vector<unsigned int> tour) {
    vector<pair<unsigned int, unsigned int>> v;

    // Pickup indices have to be less than their respective deliveries. If a pickup is found
    // with its delivery index as non-zero, that is a violation.
    vector<unsigned int> deliveries(tour.size(), 0);

    // Note that we ignore +0/-0 at the ends of the tour.
    for (unsigned int i = 1; i < tour.size() - 1; ++i)
        if (problem.has_predecessor(tour[i]))
            deliveries[tour[i]] = i;
        else {
            auto d = problem.successor_index(tour[i]);
            if (deliveries[d] > 0)
                v.push_back({deliveries[d], i});
        }

    return v;
}

void APATSPPDCallback::cut_violation(vector<unsigned int> tour, pair<unsigned int, unsigned int> index) {
    GRBLinExpr expr = 0;
    for (unsigned int i = 0; i <= index.first; ++i)
        for (unsigned int j = index.first + 1; j < tour.size(); ++j)
            expr += x[tour[i]][tour[j]] + x[tour[j]][tour[i]];

    addLazy(expr >= 3);
}
