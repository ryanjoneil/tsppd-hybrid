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
#include <tsppd/solver/oneil/oneil_tsppd_callback.h>
#include <tsppd/solver/oneil/oneil_tsppd_solver.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Solver;
using namespace std;

ONeilTSPPDSolver::ONeilTSPPDSolver(
    const TSPPDProblem& problem,
    const map<string, string> options,
    TSPSolutionWriter& writer) :
    TSPSolver(problem, options, writer),
    env(),
    model(env),
    x(),
    w(),
    start_index(problem.index("+0")),
    end_index(problem.index("-0")) {

    // Silence output.
    model.getEnv().set(GRB_IntParam_OutputFlag, 0);

    initialize_variables();
    initialize_assignment_problem_constraints();
    initialize_x_w_link_constraints();
    initialize_subtour_and_precedence_constraints();
}

TSPPDSolution ONeilTSPPDSolver::solve() {
    // Set time limit.
    if (time_limit > 0)
        model.set(GRB_DoubleParam_TimeLimit, time_limit / 1000.0);

    // Set solution limit.
    if (solution_limit > 0)
        model.set(GRB_IntParam_SolutionLimit, solution_limit);

    // Set thread count.
    model.set(GRB_IntParam_Threads, threads);

    ONeilTSPPDCallback callback(problem, x, writer);
    model.setCallback(&callback);

    model.optimize();

    // cout << endl;

    // cout << "w =" << "\t";
    // for (auto p : problem.pickup_indices())
    //     cout << setw(6) << left << problem.nodes[p] << " ";
    // cout << endl;
    // for (auto pi : problem.pickup_indices()) {
    //     cout << setw(4) << left << problem.nodes[pi] << "\t";
    //     for (auto pj : problem.pickup_indices()) {
    //         if (pi == pj)
    //             cout << "- - - ";
    //         else {
    //             for (auto w_ijk : w[{pi,pj}])
    //                 if (*model.get(GRB_DoubleAttr_X, &(w_ijk), 1) > 0.5)
    //                     cout << "1 ";
    //                 else
    //                     cout << "- ";
    //         }
    //         cout << " ";
    //     }
    //     cout << endl;
    // }
    // cout << endl;

    // cout << "x = " << endl;
    // for (unsigned int i = 0; i < problem.nodes.size(); ++i) {
    //     cout << "[" << problem.nodes[i] << "] ";
    //     for (unsigned int j = 0; j < problem.nodes.size(); ++j)
    //         if (*model.get(GRB_DoubleAttr_X, &(x[i][j]), 1) > 0.5)
    //             cout << problem.nodes[j] << " ";
    //     cout << endl;
    // }

    // cout << endl;


    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        TSPPDSolution solution(problem, get_path());
        TSPPDSearchStatistics stats(solution);
        stats.dual = model.get(GRB_DoubleAttr_ObjBound);
        writer.write(stats, true);
        return solution;
    }

    return {problem, problem.nodes};
}

void ONeilTSPPDSolver::initialize_variables() {
    for (unsigned int from = 0; from < problem.nodes.size(); ++from) {
        vector<GRBVar> x_i;

        for (unsigned int to = 0; to < problem.nodes.size(); ++to) {
            // +0 can only connect to a pickup, only deliveries can connect to -0
            auto ub = (from == to || from == end_index || to == start_index) ? 0 : 1;
            if (from == start_index && !problem.has_successor(to))
                ub = 0;
            else if (to == end_index && !problem.has_predecessor(from))
                ub = 0;
            else if (problem.has_successor(to) && problem.successor_index(to) == from)
                ub = 0;

            x_i.push_back(model.addVar(0, ub, problem.cost(from, to), GRB_BINARY));

        }

        x.push_back(x_i);
    }

    auto pickups = problem.pickup_indices();
    for (unsigned int i = 0; i < pickups.size(); ++i) {
        auto pi = pickups[i];

        for (unsigned int j = i + 1; j < pickups.size(); ++j) {
            auto pj = pickups[j];

            // auto w_ij1 = model.addVar(0, 1, 0, GRB_BINARY); // +i +j -i -j -> w_ij1
            // auto w_ij2 = model.addVar(0, 1, 0, GRB_BINARY); // +i +j -j -i -> w_ij2
            // auto w_ij3 = model.addVar(0, 1, 0, GRB_BINARY); // +i -i +j -j -> w_ij3

            // auto w_ji1 = model.addVar(0, 1, 0, GRB_BINARY); // +j +i -j -i -> w_ji1
            // auto w_ji2 = model.addVar(0, 1, 0, GRB_BINARY); // +j +i -i -j -> w_ji2
            // auto w_ji3 = model.addVar(0, 1, 0, GRB_BINARY); // +j -j +i -i -> w_ji3

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

void ONeilTSPPDSolver::initialize_assignment_problem_constraints() {
    // sum {j != i} x_ij = 1 for all i = 2,...,n
    for (unsigned int to = 0; to < problem.nodes.size(); ++to) {
        if (to == start_index)
            continue;

        GRBLinExpr expr = 0;
        for (unsigned int from = 0; from < problem.nodes.size(); ++from) {
            if (from == to)
                continue;
            expr += x[from][to];
        }

        // +0 is preceded by nothing
        auto rhs = to == start_index ? 0 : 1;
        model.addConstr(expr == rhs);
    }

    // sum {i != j} x_ij = 1 for all j = 2,...,n
    for (unsigned int from = 0; from < problem.nodes.size(); ++from) {
        if (from == start_index)
            continue;

        GRBLinExpr expr = 0;
        for (unsigned int to = 0; to < problem.nodes.size(); ++to) {
            if (from == to)
                continue;

            expr += x[from][to];
        }

        // -0 precedes nothing
        auto rhs = from == end_index ? 0 : 1;
        model.addConstr(expr == rhs);
    }
}

void ONeilTSPPDSolver::initialize_x_w_link_constraints() {
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
}

void ONeilTSPPDSolver::initialize_subtour_and_precedence_constraints() {
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

                model.addConstr(sec(pi,pj) + sec(pj,pk) + sec(pk,pi) <= 2); // +i +j +k +i
                model.addConstr(sec(pi,pj) + sec(pj,dk) + sec(dk,pi) <= 2); // +i +j -k +i
                model.addConstr(sec(pi,dj) + sec(dj,pk) + sec(pk,pi) <= 2); // +i -j +k +i
                model.addConstr(sec(pi,dj) + sec(dj,dk) + sec(dk,pi) <= 2); // +i -j -k +i

                model.addConstr(sec(di,pj) + sec(pj,pk) + sec(pk,di) <= 2); // -i +j +k -i
                model.addConstr(sec(di,pj) + sec(pj,dk) + sec(dk,di) <= 2); // -i +j -k -i
                model.addConstr(sec(di,dj) + sec(dj,pk) + sec(pk,di) <= 2); // -i -j +k -i
                model.addConstr(sec(di,dj) + sec(dj,dk) + sec(dk,di) <= 2); // -i -j -k -i

                // model.addConstr(3*(sec(pi,pj) + sec(pj,pk) + sec(pk,pi)) + x[pj][pi] + x[pk][pj] + x[pi][pk] <= 6); // +i +j +k +i
                // model.addConstr(3*(sec(pi,pj) + sec(pj,dk) + sec(dk,pi)) + x[pj][pi] + x[dk][pj] + x[pi][dk] <= 6); // +i +j -k +i
                // model.addConstr(3*(sec(pi,dj) + sec(dj,pk) + sec(pk,pi)) + x[dj][pi] + x[pk][dj] + x[pi][pk] <= 6); // +i -j +k +i
                // model.addConstr(3*(sec(pi,dj) + sec(dj,dk) + sec(dk,pi)) + x[dj][pi] + x[dk][dj] + x[pi][dk] <= 6); // +i -j -k +i

                // model.addConstr(3*(sec(di,pj) + sec(pj,pk) + sec(pk,di)) + x[pj][di] + x[pk][pj] + x[di][pk] <= 6); // -i +j +k -i
                // model.addConstr(3*(sec(di,pj) + sec(pj,dk) + sec(dk,di)) + x[pj][di] + x[dk][pj] + x[di][dk] <= 6); // -i +j -k -i
                // model.addConstr(3*(sec(di,dj) + sec(dj,pk) + sec(pk,di)) + x[dj][di] + x[pk][dj] + x[di][pk] <= 6); // -i -j +k -i
                // model.addConstr(3*(sec(di,dj) + sec(dj,dk) + sec(dk,di)) + x[dj][di] + x[dk][dk] + x[di][dk] <= 6); // -i -j -k -i
            }
        }
    }
}

GRBLinExpr ONeilTSPPDSolver::sec(unsigned int i, unsigned int j) {
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

vector<unsigned int> ONeilTSPPDSolver::get_path() {
    vector<unsigned int> path{0};
    unsigned int from = 0;
    for (unsigned int i = 0; i < problem.nodes.size() - 1; ++i) {
        for (unsigned int to = 0; to < problem.nodes.size(); ++to) {
            if (*model.get(GRB_DoubleAttr_X, &(x[from][to]), 1) > 0.5) {
                path.push_back(to);
                from = to;
                break;
            }
        }
    }
    return path;
}
