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

#include <tsppd/solver/focacci/focacci_tsppd_solver.h>
#include <tsppd/solver/focacci/focacci_tsppd_space.h>
#include <tsppd/util/exception.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Solver;
using namespace TSPPD::Util;
using namespace std;

FocacciTSPPDSolver::FocacciTSPPDSolver(
    const TSPPDProblem& problem,
    const map<string, string> options,
    TSPSolutionWriter& writer) :
    FocacciTSPSolver(problem, options, writer) {

    initialize_tsppd_options();
}

void FocacciTSPPDSolver::initialize_tsppd_options() {
    // Precedence propagation
    if (options["precede"] == "" || options["precede"] == "set")
        precede_type = PRECEDE_SET;
    else if (options["precede"] == "cost")
        precede_type = PRECEDE_COST;
    else if (options["precede"] == "all")
        precede_type = PRECEDE_ALL;
    else
        throw TSPPDException("precede can be either set, cost, or all");

    // Additive bounding reduced cost fixing
    filter_add = false;
    filter_ap = false;
    if (options["filter"] == "add")
        filter_add = true;
    else if (options["filter"] == "ap")
        filter_ap = true;
    else if (options["filter"] != "" && options["filter"] != "none")
        throw TSPPDException("filter can be either add, ap, or none");

    // Order Matching Constraint propagation
    if (options["omc"] == "" || options["omc"] == "off")
        omc = false;
    else if (options["omc"] == "on")
        omc = true;
    else
        throw TSPPDException("omc can be either on or off");
}

shared_ptr<FocacciTSPSpace> FocacciTSPPDSolver::build_space() {
    auto space = make_shared<FocacciTSPPDSpace>(problem);
    space->initialize_precedence_propagators(precede_type);
    if (filter_add)
        space->initialize_additive_bounding();
    else if (filter_ap)
        space->initialize_assignment_propagator();
    if (omc)
        space->initialize_omc_constraints();
    return space;
}
