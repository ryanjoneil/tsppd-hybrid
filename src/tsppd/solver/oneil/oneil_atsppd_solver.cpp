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

#include <cstdio>
#include <iomanip>
#include <vector>

#include <tsppd/data/tsppd_search_statistics.h>
#include <tsppd/solver/ap/ap_atsppd_callback.h>
#include <tsppd/solver/oneil/oneil_atsppd_solver.h>
#include <tsppd/util/exception.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Solver;
using namespace TSPPD::Util;
using namespace std;

ONeilATSPPDSolver::ONeilATSPPDSolver(
    const TSPPDProblem& problem,
    const map<string, string> options,
    TSPSolutionWriter& writer) :
    APATSPPDSolver(problem, options, writer),
    w() {

    initialize_oneil_options();
    initialize_oneil_variables();
    initialize_oneil_constraints();
}

TSPPDSolution ONeilATSPPDSolver::solve() {
    configure_solver();

    APATSPPDCallback callback(problem, x, sec, writer);
    model.setCallback(&callback);

    model.optimize();
    return solution();
}

void ONeilATSPPDSolver::initialize_oneil_options() {
    if (options["relax"] == "" || options["relax"] == "off")
        relaxed = false;
    else if (options["relax"] == "on")
        relaxed = true;
    else
        throw TSPPDException("relax can be either on or off");

    if (options["sec"] == "" || options["sec"] == "subtour")
        sec = ATSP_SEC_SUBTOUR;
    else if (options["sec"] == "cutset")
        sec = ATSP_SEC_CUTSET;
    else
        throw TSPPDException("sec can be either subtour or cutset");
}

void ONeilATSPPDSolver::initialize_oneil_variables() {
    auto pickups = problem.pickup_indices();
    for (unsigned int i = 0; i < pickups.size(); ++i) {
        auto pi = pickups[i];

        for (unsigned int j = i + 1; j < pickups.size(); ++j) {
            auto pj = pickups[j];

            auto w_ij1 = model.addVar(0, 1, 0, GRB_CONTINUOUS); // +i +j -i -j -> w_ij1
            auto w_ij2 = model.addVar(0, 1, 0, GRB_CONTINUOUS); // +i +j -j -i -> w_ij2
            auto w_ij3 = model.addVar(0, 1, 0, GRB_CONTINUOUS); // +i -i +j -j -> w_ij3

            auto w_ji1 = model.addVar(0, 1, 0, GRB_CONTINUOUS); // +j +i -j -i -> w_ji1
            auto w_ji2 = model.addVar(0, 1, 0, GRB_CONTINUOUS); // +j +i -i -j -> w_ji2
            auto w_ji3 = model.addVar(0, 1, 0, GRB_CONTINUOUS); // +j -j +i -i -> w_ji3

            w[{pi, pj}] = {w_ij1, w_ij2, w_ij3};
            w[{pj, pi}] = {w_ji1, w_ji2, w_ji3};

            model.addConstr(w_ij1 + w_ij2 + w_ij3 + w_ji1 + w_ji2 + w_ji3 == 1);
        }
    }
}

void ONeilATSPPDSolver::initialize_oneil_constraints() {
    // Link x and w variables
    for (unsigned int i = 0; i < x.size(); ++i) {
        if (i == start_index || i == end_index)
            continue;

        for (unsigned int j = 0; j < x.size(); ++j) {
            if (j == start_index || j == end_index || j == i)
                continue;

            auto x_ij = x[i][j];

            if (problem.has_successor(i) && problem.has_successor(j)) {
                // x(+i,+j) <= w_ij1 + w_ij2
                auto pi = i;
                auto pj = j;
                model.addConstr(x_ij <= w[{pi,pj}][0] + w[{pi,pj}][1]);

                // x(+i,-i) <= w_ij3 + w_ji2 + w_ji3
                auto x_ii = x[pi][problem.successor_index(i)];
                model.addConstr(x_ii <= w[{pi,pj}][2] + w[{pj,pi}][1] + w[{pj,pi}][2]);

            } else if (problem.has_successor(i) && problem.has_predecessor(j)) {
                // x(+i,-j) <= w_ji1
                auto pi = i;
                auto pj = problem.predecessor_index(j);
                if (pi != pj)
                    model.addConstr(x_ij <= w[{pj,pi}][0]);

            } else if (problem.has_predecessor(i) && problem.has_successor(j)) {
                // x(-i,+j) <= w_ij3
                auto pi = problem.predecessor_index(i);
                auto pj = j;
                if (pi != pj)
                    model.addConstr(x_ij <= w[{pi,pj}][2]);

            } else if (problem.has_predecessor(i) && problem.has_predecessor(j)) {
                // x(-i,-j) <= w_ij1 + w_ji2
                auto pi = problem.predecessor_index(i);
                auto pj = problem.predecessor_index(j);
                model.addConstr(x_ij <= w[{pi,pj}][0] + w[{pj,pi}][1]);
            }
        }
    }

    // Subtour elimination and precedence constraints
    if (!relaxed) {
        for (auto pi : problem.pickup_indices()) {
            auto di = problem.successor_index(pi);

            for (auto pj : problem.pickup_indices()) {
                if (pi == pj)
                    continue;

                auto dj = problem.successor_index(pj);
                for (auto pk : problem.pickup_indices()) {
                    if (pk == pi || pk == pj)
                        continue;

                    auto dk = problem.successor_index(pk);

                    model.addConstr(s(pi,pj) + s(pj,pk) + s(pk,pi) <= 2); // +i +j +k +i
                    model.addConstr(s(pi,pj) + s(pj,dk) + s(dk,pi) <= 2); // +i +j -k +i
                    model.addConstr(s(pi,dj) + s(dj,pk) + s(pk,pi) <= 2); // +i -j +k +i
                    model.addConstr(s(pi,dj) + s(dj,dk) + s(dk,pi) <= 2); // +i -j -k +i

                    model.addConstr(s(di,pj) + s(pj,pk) + s(pk,di) <= 2); // -i +j +k -i
                    model.addConstr(s(di,pj) + s(pj,dk) + s(dk,di) <= 2); // -i +j -k -i
                    model.addConstr(s(di,dj) + s(dj,pk) + s(pk,di) <= 2); // -i -j +k -i
                    model.addConstr(s(di,dj) + s(dj,dk) + s(dk,di) <= 2); // -i -j -k -i
                }
            }
        }
    }
}

GRBLinExpr ONeilATSPPDSolver::s(unsigned int i, unsigned int j) {
    GRBLinExpr expr = 0;

    if (problem.has_successor(i) && problem.has_successor(j)) {
        auto pi = i;
        auto pj = j;
        expr = w[{pi,pj}][0] + w[{pi,pj}][1] + w[{pi,pj}][2];

    } else if (problem.has_successor(i) && problem.has_predecessor(j)) {
        auto pi = i;
        auto pj = problem.predecessor_index(j);
            if (pi != pj)
                expr = w[{pi,pj}][0] + w[{pi,pj}][1] + w[{pi,pj}][2] + w[{pj,pi}][0] + w[{pj,pi}][1];

    } else if (problem.has_predecessor(i) && problem.has_successor(j)) {
        auto pi = problem.predecessor_index(i);
        auto pj = j;
        if (pi != pj)
            expr = w[{pi,pj}][2];

    } else if (problem.has_predecessor(i) && problem.has_predecessor(j)) {
        auto pi = problem.predecessor_index(i);
        auto pj = problem.predecessor_index(j);
        if (pi != pj)
            expr = w[{pi,pj}][0] + w[{pi,pj}][2] + w[{pj,pi}][1];
    }

    return expr;
}
