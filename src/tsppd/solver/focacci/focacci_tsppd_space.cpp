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

#include <tsppd/solver/focacci/propagator/focacci_tsppd_additive_propagator.h>
#include <tsppd/solver/focacci/propagator/focacci_tsppd_assignment_propagator.h>
#include <tsppd/solver/focacci/propagator/focacci_tsppd_heldkarp_propagator.h>
#include <tsppd/solver/focacci/propagator/focacci_tsppd_omc_propagator.h>
#include <tsppd/solver/focacci/propagator/focacci_tsppd_precede_cost_propagator.h>
#include <tsppd/solver/focacci/propagator/focacci_tsppd_precede_set_propagator.h>
#include <tsppd/solver/focacci/focacci_tsppd_space.h>

using namespace Gecode;
using namespace TSPPD::Data;
using namespace TSPPD::Solver;
using namespace std;

FocacciTSPPDSpace::FocacciTSPPDSpace(const TSPPDProblem& problem) :
    FocacciTSPSpace(problem) { }

FocacciTSPPDSpace::FocacciTSPPDSpace(FocacciTSPPDSpace& s) :
    FocacciTSPSpace(s) { }

Gecode::Space* FocacciTSPPDSpace::copy() {
    return new FocacciTSPPDSpace(*this);
}

void FocacciTSPPDSpace::initialize_constraints() {
    FocacciTSPSpace::initialize_constraints();

    // Constraints: Always start with +0 and end with -0.
    unsigned int start_index = 0;
    unsigned int end_index = problem.index(problem.successor(problem.nodes[0]));
    rel(*this, next[end_index] == start_index);

    // Must start with a pickup and end with a delivery.
    for (auto node : problem.nodes) {
        auto i = problem.index(node);

        if (i == start_index || i == end_index)
            continue;
        else if (problem.has_successor(node))
            // Pickup cannot directly precede the end of the route.
            rel(*this, next[i] != end_index);
        else if (problem.has_predecessor(node))
            // Delivery cannot directly succeed the start of the route.
            rel(*this, next[start_index] != i);
    }

    // Deliveries cannot precede pickups.
    for (auto node : problem.nodes) {
        if (problem.has_successor(node)) {
            auto pickup = problem.index(node);
            auto delivery = problem.index(problem.successor(node));
            if (pickup != start_index)
                rel(*this, next[delivery] != pickup);
        }
    }
}

void FocacciTSPPDSpace::initialize_precedence_propagators(const FocacciTSPPDPrecedePropagatorType precede_type) {
    if (precede_type == PRECEDE_SET || precede_type == PRECEDE_ALL)
        tsppd_precede_set(*this, next, problem);
    if (precede_type == PRECEDE_COST || precede_type == PRECEDE_ALL)
        tsppd_precede_cost(*this, length, next, problem);
}

void FocacciTSPPDSpace::initialize_additive_propagator() {
    tsppd_additive(*this, next, length, problem);
}

void FocacciTSPPDSpace::initialize_assignment_propagator() {
    cout << "foo\n";
    tsppd_assignment(*this, next, length, problem);
}

void FocacciTSPPDSpace::initialize_heldkarp_propagator() {
    tsppd_heldkarp(*this, next, length, problem);
}

void FocacciTSPPDSpace::initialize_omc_constraints() {
    tsppd_omc(*this, next, problem);
}
