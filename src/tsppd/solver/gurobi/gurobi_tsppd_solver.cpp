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

#include <tsppd/solver/gurobi/gurobi_tsppd_solver.h>
#include <tsppd/solver/gurobi/callback/gurobi_precedence_callback.h>
#include <tsppd/util/exception.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Solver;
using namespace TSPPD::Util;
using namespace std;

GurobiTSPPDSolver::GurobiTSPPDSolver(
    const TSPPDProblem& problem,
    const std::map<std::string, std::string> options,
    TSPSolutionWriter& writer) :
    GurobiTSPSolver(problem, options, writer) {

    initialize_tsppd_options();
    initialize_tsppd_constraints();
    if (omc)
        initialize_omc_constraints();
}

void GurobiTSPPDSolver::initialize_tsppd_options() {
    if (options["omc"] == "" || options["omc"] == "off")
        omc = false;
    else if (options["omc"] == "on")
        omc = true;
    else
        throw TSPPDException("omc can be either on or off");

    if (options["omc-lazy"] == "" || options["omc-lazy"] == "off")
        omc_lazy = false;
    else if (options["omc-lazy"] == "on")
    omc_lazy = true;
    else
        throw TSPPDException("omc-lazy can be either on or off");
}

void GurobiTSPPDSolver::initialize_tsppd_constraints() {
    // +0 and -0 are connected.
    unsigned int start = 0;
    auto end = problem.successor_index(start);
    if (problem.has_successor(start))
        model.addConstr(arcs[{start, end}] == 1);
    else
        return;

    // Must start with a pickup and end with a delivery.
    GRBLinExpr expr_pickup_start = 0;
    GRBLinExpr expr_delivery_end = 0;
    for (unsigned int node = 0; node < problem.nodes.size(); ++node) {
        if (node == start || node == end)
            continue;

        if (problem.has_successor(node))
            // Is eligible for start of the tour.
            expr_pickup_start += arcs[{start, node}];
        else if (problem.has_predecessor(node))
            // Is eligible for end of the tour.
            expr_delivery_end += arcs[{node, end}];
    }

    model.addConstr(expr_pickup_start == 1);
    model.addConstr(expr_delivery_end == 1);
}

void GurobiTSPPDSolver::initialize_omc_constraints() {
    auto pickups = problem.pickup_indices();

    for (size_t p1_idx = 0; p1_idx < pickups.size(); ++p1_idx) {
        auto p1 = pickups[p1_idx];
        auto d1 = problem.successor_index(p1);

        for (size_t p2_idx = p1_idx + 1; p2_idx < pickups.size(); ++p2_idx) {
            auto p2 = pickups[p2_idx];
            auto d2 = problem.successor_index(p2);

            // (+i -i) + (+j -j) + (+i +j)
            model.addConstr(arcs[{p1, d1}] + arcs[{p2, d2}] + arcs[{p1, p2}] <= 2);

            // (+i -i) + (+j -j) + (-i -j)
            model.addConstr(arcs[{p1, d1}] + arcs[{p2, d2}] + arcs[{d1, d2}] <= 2);

            // (+i -j) + (+j -i) + (+i +j)
            model.addConstr(arcs[{p1, d2}] + arcs[{p2, d1}] + arcs[{p1, p2}] <= 2);

            // (+i -j) + (+j -i) + (-i ij)
            model.addConstr(arcs[{p1, d2}] + arcs[{p2, d1}] + arcs[{d1, d2}] <= 2);
        }
    }
}

void GurobiTSPPDSolver::initialize_callbacks() {
    callbacks.push_back(
        static_cast<shared_ptr<GurobiTSPCallback>>(
            make_shared<GurobiPrecedenceCallback>(
                options["sec"],
                problem,
                arcs,
                omc_lazy
            )
        )
    );
}

