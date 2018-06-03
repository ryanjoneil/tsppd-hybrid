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
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include <boost/algorithm/string/join.hpp>

#include <tsppd/io/tsp_solution_writer.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace std;

TSPSolutionWriter::TSPSolutionWriter(
    const TSPPDProblem& problem,
    const string solver,
    const map<string, string> options,
    const TSPSolutionFormat format) :
    problem(problem),
    solver(solver),
    options(options),
    format(format),
    start(clock()) { }

void TSPSolutionWriter::write_header() {
    if (format == HUMAN) {
        cout << "instance         size   solver        elapsed   dual      primal    nodes     fails     depth     ";
        for (auto opt : options)
            cout << setfill(' ') << setw(10) << left << opt.first;
        cout << endl;

        for (unsigned int i = 0; i < TSPWriterSeparatorLength; ++i)
            cout << '=';
        cout << endl;

    } else if (format == CSV) {
        cout << "instance,size,solver,elapsed,dual,primal,nodes,fails,depth";
        for (auto opt : options)
            cout << "," << opt.first;
        cout << "," << "tour" << endl;
    }
}

void TSPSolutionWriter::write(const TSPPDSearchStatistics& stats, const bool force) {
    auto now = clock();
    auto elapsed = ((double) (now - start)) / CLOCKS_PER_SEC;

    auto dual_str = (stats.has_dual() ? to_string(stats.dual) : "");
    auto primal_str = (stats.has_primal() ? to_string(stats.primal) : "");
    auto elapsed_str = to_string(elapsed);

    if (format == HUMAN) {
        if (!force && last_dual_str == dual_str && last_primal_str == primal_str)
            return;

        last_dual_str = dual_str;
        last_primal_str = primal_str;

        stringstream s;
        s << fixed << setprecision(3) << elapsed;
        elapsed_str = s.str();
    }

    std::vector<std::string> row {
        problem.name,
        to_string(problem.nodes.size()),
        solver,
        elapsed_str,
        dual_str,
        primal_str,
        stats.has_nodes() ? to_string(stats.nodes) : "",
        stats.has_fails() ? to_string(stats.fails) : "",
        stats.has_depth() ? to_string(stats.depth) : ""
    };

    for (auto opt : options)
        row.push_back(opt.second);

    if (format == HUMAN) {
        cout << setfill(' ') << setw(17) << left << row[0]
             << setfill(' ') << setw(7) << left << row[1]
             << setfill(' ') << setw(14) << left << row[2];
        for (unsigned int i = 3; i < row.size(); ++i)
            cout << setfill(' ') << setw(10) << left << row[i];
        cout << endl;
    } else if (format == CSV) {
        row.push_back(boost::algorithm::join(stats.tour, " "));
        cout << boost::algorithm::join(row, ",") << endl;
    }
}
