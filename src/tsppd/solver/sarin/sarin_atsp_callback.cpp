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
