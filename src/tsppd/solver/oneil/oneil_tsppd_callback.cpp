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

#include <tsppd/data/tsppd_search_statistics.h>
#include <tsppd/data/tsppd_solution.h>
#include <tsppd/solver/oneil/oneil_tsppd_callback.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Solver;
using namespace std;

ONeilTSPPDCallback::ONeilTSPPDCallback(
    const TSPPDProblem& problem,
    vector<vector<GRBVar>> x,
    map<pair<unsigned int, unsigned int>, vector<GRBVar>> w,
    TSPSolutionWriter& writer) :
    problem(problem), x(x), w(w), writer(writer) { }

void ONeilTSPPDCallback::callback() {
    if (where == GRB_CB_MIP) {
        // Just log bounds.
        TSPPDSearchStatistics stats;
        stats.primal = getDoubleInfo(GRB_CB_MIP_OBJBST);
        stats.dual = max(0.0, getDoubleInfo(GRB_CB_MIP_OBJBND));
        writer.write(stats);

    } else if (where == GRB_CB_MIPSOL) {
        auto s = subtours();

        if (s.size() > 1) {
            // Elinate subtours.
            for (auto subtour : s)
                cut_subtour(subtour);

        } else {
            // We have a single tour. Check it for precedence violations.
            auto v = violations(s.front());
            if (v.size() > 0)
                for (auto vi : v)
                    cut_violation(s.front(), vi);

            else {
                // Single tour with no precedence violations found.
                TSPPDSolution solution(problem, s.front());
                TSPPDSearchStatistics stats(solution);
                stats.dual = max(0.0, getDoubleInfo(GRB_CB_MIPSOL_OBJBND));
                writer.write(stats);
            }
        }
    }
}

vector<vector<unsigned int>> ONeilTSPPDCallback::subtours() {
    vector<vector<unsigned int>> s;

    // Set of unseen nodes.
    set<unsigned int> unseen;
    for (unsigned int i = 0; i < problem.nodes.size(); ++i)
        unseen.insert(i);

    while (!unseen.empty()) {
        unsigned int next = *unseen.begin();
        unseen.erase(next);

        vector<unsigned int> subtour{next};

        auto done = false;
        while (!done) {
            done = true;
            for (unsigned int i = 0; i < problem.nodes.size(); ++i) {
                if (getSolution(x[next][i]) > 0.5) {
                    subtour.push_back(i);
                    if (unseen.find(i) != unseen.end()) {
                        done = false;
                        unseen.erase(i);
                        next = i;
                    }
                    break;
                }
            }
        }

        s.push_back(subtour);
    }

    return s;
}

void ONeilTSPPDCallback::cut_subtour(const vector<unsigned int>& subtour) {
    set<unsigned int> S(subtour.begin(), subtour.end());

    vector<unsigned int> T;
    for (unsigned int i = 0; i < problem.nodes.size(); ++i)
        if (S.find(i) == S.end())
            T.push_back(i);

    GRBLinExpr expr = 0;
    for (auto n1 : S)
        for (auto n2 : T)
            expr += x[n1][n2] + x[n2][n1];

    addLazy(expr >= 1);
}


vector<pair<unsigned int, unsigned int>>ONeilTSPPDCallback::violations(vector<unsigned int> tour) {
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

void ONeilTSPPDCallback::cut_violation(vector<unsigned int> tour, pair<unsigned int, unsigned int> index) {
    GRBLinExpr expr = 0;
    for (unsigned int i = 0; i <= index.first; ++i)
        for (unsigned int j = index.first + 1; j < tour.size(); ++j)
            expr += x[tour[i]][tour[j]] + x[tour[j]][tour[i]];

    addLazy(expr >= 3);
}
