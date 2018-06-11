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

#include <iomanip>
#include <set>

#include <tsppd/data/tsppd_search_statistics.h>
#include <tsppd/data/tsppd_solution.h>
#include <tsppd/solver/sarin/sarin_tsp_callback.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Solver;
using namespace std;

SarinTSPCallback::SarinTSPCallback(
    const TSPPDProblem& problem,
    vector<vector<GRBVar>> x,
    vector<vector<GRBVar>> y,
    TSPSolutionWriter& writer) :
    problem(problem), x(x), y(y), writer(writer) { }

void SarinTSPCallback::callback() {
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

vector<vector<unsigned int>> SarinTSPCallback::subtours() {
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

void SarinTSPCallback::cut_subtour(vector<unsigned int> subtour) {
    if (subtour.front() != subtour.back())
        return;

    cout << "subtour: " << endl;
    for (auto i : subtour)
        cout << problem.nodes[i] << " ";
    cout << endl;

    auto i = 0; // start and end of the subtour
    for (unsigned int j = 1; j < subtour.size() - 2; ++j)
        for (unsigned int k = j + 1; k < subtour.size() - 1; ++k)
            addLazy(
                y[subtour[i]][subtour[j]] +
                x[subtour[j]][subtour[i]] +
                y[subtour[j]][subtour[k]] +
                y[subtour[k]][subtour[i]] <= 2
            );

    cout << endl;
}


vector<pair<unsigned int, unsigned int>>SarinTSPCallback::violations(vector<unsigned int> tour) {
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

void SarinTSPCallback::cut_violation(vector<unsigned int> tour, pair<unsigned int, unsigned int> index) {
    // This method takes in partial tours of the form:
    //
    //     -i ... j ... k ... l ... +i
    //
    // It adds cuts of the form:
    //
    //    y(-i,j) + x(j,-i) + y(j,k) + y(k,+i) <= 2
    //    etc.
    auto d = index.first;
    auto p = index.second;

    cout << "cut violation: " << problem.nodes[tour[d]] << " < " << problem.nodes[tour[p]] << endl;

    GRBLinExpr expr = y[tour[p]][tour[d]];
    for (unsigned int j = d; j < p; ++j)
        expr += y[tour[j]][tour[j+1]];
    addLazy(expr <= index.second - index.first);

}
