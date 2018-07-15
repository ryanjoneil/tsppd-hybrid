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
#include <tsppd/solver/ap/ap_atsp_callback.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Solver;
using namespace std;

APATSPCallback::APATSPCallback(
    const TSPPDProblem& problem,
    vector<vector<GRBVar>> x,
    const ATSPSECType sec_type,
    TSPSolutionWriter& writer) :
    problem(problem), x(x), sec_type(sec_type), writer(writer), primal(-1) { }

void APATSPCallback::callback() {
    if (where == GRB_CB_MIP)
        log_mip();

    else if (where == GRB_CB_MIPSOL) {
        auto s = subtours();
        if (s.size() > 1)
            for (auto subtour : s)
                cut_subtour(subtour);
        else
            log_mipsol(s.front());
    }
}

void APATSPCallback::log_mip() {
    TSPPDSearchStatistics stats;

    // For some reason this moodel produces invalid primal bounds during search sometimes,
    // so we instead get primal bounds directly off of MIPSOL callbacks.
    if (primal >= 0)
        stats.primal = primal;
    else
        stats.primal = getDoubleInfo(GRB_CB_MIP_OBJBST);

    stats.dual = max(0.0, getDoubleInfo(GRB_CB_MIP_OBJBND));
    writer.write(stats);
}

void APATSPCallback::log_mipsol(vector<unsigned int>& tour) {
    TSPPDSolution solution(problem, tour);
    TSPPDSearchStatistics stats(solution);
    primal = solution.cost;
    stats.dual = max(0.0, getDoubleInfo(GRB_CB_MIPSOL_OBJBND));
    writer.write(stats);
}

vector<vector<unsigned int>> APATSPCallback::subtours() {
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
                    if (unseen.find(i) != unseen.end()) {
                        subtour.push_back(i);
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

void APATSPCallback::cut_subtour(const vector<unsigned int>& subtour) {
    if (sec_type == ATSP_SEC_CUTSET)
        cut_subtour_cutset(subtour);
    else if (sec_type == ATSP_SEC_SUBTOUR)
        cut_subtour_subtour(subtour);
}

void APATSPCallback::cut_subtour_cutset(const vector<unsigned int>& subtour) {
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

void APATSPCallback::cut_subtour_subtour(const vector<unsigned int>& subtour) {
    GRBLinExpr expr = 0;

    for (unsigned int i = 0; i < subtour.size(); ++i)
        for (unsigned int j = 0 ; j < subtour.size(); ++j)
            if (i != j)
                expr += x[subtour[i]][subtour[j]];

    addLazy(expr <= subtour.size() - 1);
}
