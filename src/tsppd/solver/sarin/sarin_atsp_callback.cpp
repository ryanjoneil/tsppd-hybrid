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
#include <tsppd/solver/sarin/sarin_atsp_callback.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Solver;
using namespace std;

SarinATSPCallback::SarinATSPCallback(
    const TSPPDProblem& problem,
    vector<vector<GRBVar>> x,
    map<pair<unsigned int, unsigned int>, GRBVar> y,
    const ATSPSECType sec,
    TSPSolutionWriter& writer) :
    APATSPCallback(problem, x, sec, writer),
    y(y),
    start_index(problem.index("+0")),
    end_index(problem.index("-0")) { }

void SarinATSPCallback::callback() {
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

void SarinATSPCallback::cut_subtour(const vector<unsigned int>& subtour) {
    if (sec_type == ATSP_SEC_OTHER)
        cut_subtour_other(subtour);
    else
        APATSPCallback::cut_subtour(subtour);
}

void SarinATSPCallback::cut_subtour_other(const vector<unsigned int>& subtour) {
    if (subtour.front() == start_index && subtour.back() == end_index) {
        // +0 ... i ... j ... -0
        return;

    } else {
        // i ... j ... k ... i
        auto i = 0; // start and end of the subtour
        for (unsigned int j = 1; j < subtour.size() - 1; ++j) {
            for (unsigned int k = j + 1; k < subtour.size(); ++k) {
                auto ni = subtour[i];
                auto nj = subtour[j];
                auto nk = subtour[k];
                addLazy(y[{ni,nj}] + x[nj][ni] + y[{nj,nk}] + y[{nk,ni}] <= 2);
            }
        }
    }
}


// vector<pair<unsigned int, unsigned int>>SarinTSPCallback::violations(vector<unsigned int> tour) {
//     vector<pair<unsigned int, unsigned int>> v;

//     // Pickup indices have to be less than their respective deliveries. If a pickup is found
//     // with its delivery index as non-zero, that is a violation.
//     vector<unsigned int> deliveries(tour.size(), 0);

//     // Note that we ignore +0/-0 at the ends of the tour.
//     for (unsigned int i = 1; i < tour.size() - 1; ++i)
//         if (problem.has_predecessor(tour[i]))
//             deliveries[tour[i]] = i;
//         else {
//             auto d = problem.successor_index(tour[i]);
//             if (deliveries[d] > 0)
//                 v.push_back({deliveries[d], i});
//         }

//     return v;
// }

// void SarinTSPCallback::cut_violation(vector<unsigned int> tour, pair<unsigned int, unsigned int> index) {
//     // This method takes in partial tours of the form:
//     //
//     //     -i ... j ... k ... l ... +i
//     //
//     // It adds cuts of the form:
//     //
//     //    y(-i,j) + x(j,-i) + y(j,k) + y(k,+i) <= 2
//     //    etc.
//     auto d = index.first;
//     auto p = index.second;

//     // cout << "cut violation: " << problem.nodes[tour[d]] << " < " << problem.nodes[tour[p]] << endl;

//     GRBLinExpr expr = y[tour[p]][tour[d]];
//     for (unsigned int j = d; j < p; ++j)
//         expr += y[tour[j]][tour[j+1]];
//     addLazy(expr <= index.second - index.first);

// }
