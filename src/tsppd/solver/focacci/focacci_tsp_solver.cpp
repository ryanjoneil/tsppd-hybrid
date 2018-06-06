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

#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <gecode/gist.hh>
#include <gecode/search.hh>

#include <tsppd/data/tsppd_search_statistics.h>
#include <tsppd/solver/focacci/focacci_tsp_solver.h>
#include <tsppd/solver/focacci/focacci_tsp_space.h>
#include <tsppd/util/exception.h>

using namespace Gecode;
using namespace Gecode::Search;
using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Solver;
using namespace TSPPD::Util;
using namespace std;

FocacciTSPSolver::FocacciTSPSolver(
    const TSPPDProblem& problem,
    const map<string, string> options,
    TSPSolutionWriter& writer) :
    TSPSolver(problem, options, writer),
    discrepancy_limit(0) {

    initialize_tsp_options();
}

TSPPDSolution FocacciTSPSolver::solve() {
    auto space = build_space();
    space->initialize_constraints();
    space->initialize_dual(dual_type);
    space->initialize_brancher(brancher_type);

    vector<string> best_order = problem.nodes;
    auto best_cost = numeric_limits<int>::max();

#ifdef GIST
    if (gist) {
        Gist::Options o;
        Gist::Print<FocacciTSPSpace> p("next");
        o.inspect.click(&p);

        if (search_engine == SEARCH_DFS)
            Gist::dfs(space.get(), o);
        else
            Gist::bab(space.get(), o);
        return TSPPDSolution(problem, problem.nodes);
    }
#endif

    Options o;
    o.threads = threads;

    if (time_limit > 0)
        o.stop = new TimeStop(time_limit);

    // Discrepacy limit
    if (discrepancy_limit > 0)
        o.d_l = discrepancy_limit;

    unique_ptr<Base<FocacciTSPSpace>> engine;
    if (search_engine == SEARCH_DFS)
        engine = make_unique<DFS<FocacciTSPSpace>>(space.get(), o);
    else if (search_engine == SEARCH_LDS)
        engine = make_unique<LDS<FocacciTSPSpace>>(space.get(), o);
    else
        engine = make_unique<BAB<FocacciTSPSpace>>(space.get(), o);

    bool has_solution_limit = (solution_limit > 0);
    while (auto s = unique_ptr<FocacciTSPSpace>(engine->next())) {
        auto order = s->solution();
        auto cost = s->cost().val();

        TSPPDSolution solution(problem, order);

        auto gecode_stats = engine->statistics();
        TSPPDSearchStatistics stats(solution);
        stats.nodes = gecode_stats.node;
        stats.fails = gecode_stats.fail;
        stats.depth = gecode_stats.depth;

        writer.write(stats);

        if (cost < best_cost) {
            best_cost = cost;
            best_order = order;
        }

        if (has_solution_limit && --solution_limit <= 0)
            break;
    }

    TSPPDSolution solution(problem, best_order);
    auto gecode_stats = engine->statistics();
    TSPPDSearchStatistics stats(solution);
    stats.nodes = gecode_stats.node;
    stats.fails = gecode_stats.fail;
    stats.depth = gecode_stats.depth;

    writer.write(stats, true);
    return solution;
}

void FocacciTSPSolver::initialize_tsp_options() {
    initialize_option_brancher();
    initialize_option_discrepancy_limit();
    initialize_option_dual_bound();
    initialize_option_gist();
    initialize_option_search();
}

void FocacciTSPSolver::initialize_option_brancher() {
    brancher_type = BRANCHER_REGRET;
    auto brancher_pair = options.find("brancher");
    if (brancher_pair != options.end()) {
        if (brancher_pair->second == "cn")
            brancher_type = BRANCHER_CN;
        else if (brancher_pair->second == "seq-cn")
            brancher_type = BRANCHER_SEQ_CN;
        else if (brancher_pair->second != "regret")
            throw TSPPDException("invalid brancher type '" + brancher_pair->second + "'");
    }
}

void FocacciTSPSolver::initialize_option_discrepancy_limit() {
    discrepancy_limit = 0;
    auto discrepancy_limit_pair = options.find("dl");
    if (discrepancy_limit_pair != options.end()) {
        try {
            discrepancy_limit = stoi(discrepancy_limit_pair->second);
         } catch (exception &e) {
            throw TSPPDException("discrepancy limit must be an integer");
         }
        if (discrepancy_limit < 1)
            throw TSPPDException("discrepancy limit must be >= 1");
    }
}

void FocacciTSPSolver::initialize_option_dual_bound() {
    dual_type = DUAL_NONE;
    auto dual_pair = options.find("dual");
    if (dual_pair != options.end()) {
        if (dual_pair->second == "cn")
            dual_type = DUAL_CN;
        else if (dual_pair->second != "none")
            throw TSPPDException("invalid dual type '" + dual_pair->second + "'");
    }
}

void FocacciTSPSolver::initialize_option_gist() {
    gist = options.find("gist") != options.end();
}

void FocacciTSPSolver::initialize_option_search() {
    search_engine = SEARCH_BAB;
    auto search_pair = options.find("search");
    if (search_pair != options.end()) {
        if (search_pair->second == "dfs")
            search_engine = SEARCH_DFS;
        else if (search_pair->second == "lds")
            search_engine = SEARCH_LDS;
        else if (search_pair->second != "bab")
            throw TSPPDException("invalid search engine '" + search_pair->second + "'");
    }
}

shared_ptr<FocacciTSPSpace> FocacciTSPSolver::build_space() {
    return make_shared<FocacciTSPSpace>(problem);
}
