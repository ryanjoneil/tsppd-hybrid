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

#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include <tsppd/data/tsppd_problem.h>
#include <tsppd/data/tsppd_problem_generator.h>
#include <tsppd/io/tsp_problem_reader.h>
#include <tsppd/io/tsp_problem_writer.h>
#include <tsppd/solver/enumerative/enumerative_tsp_solver.h>
#include <tsppd/solver/enumerative/enumerative_tsppd_solver.h>
#include <tsppd/solver/focacci/focacci_tsp_solver.h>
#include <tsppd/solver/focacci/focacci_tsppd_solver.h>
#include <tsppd/solver/ruland/ruland_tsp_solver.h>
#include <tsppd/solver/ruland/ruland_tsppd_plus_solver.h>
#include <tsppd/solver/ruland/ruland_tsppd_solver.h>
#include <tsppd/solver/sarin/sarin_tsp_solver.h>
#include <tsppd/solver/sarin/sarin_tsppd_solver.h>
#include <tsppd/solver/tsp_solver.h>
#include <tsppd/util/exception.h>
#include <tsppd/util/stacktrace.h>

using namespace TSPPD::Solver;
using namespace std;
namespace po = boost::program_options;

int main(int argc, char** argv) {
    TSPPD::Util::enable_stacktraces();

    po::options_description desc("tsppd command line interface");
    desc.add_options()
        ("help,h", "produce help message")
        ("no-header,H", "do not print csv header line")
        ("solver,s", po::value<string>(), "{tsp|tsppd}-{enum|focacci|ruland|sarin}")
        ("input,i", po::value<string>(), "input tsplib file")
        ("format,f", po::value<string>(), "output format: {human|csv} (default=human)")
        ("random-seed,r", po::value<unsigned int>(), "random seed (default=0)")
        ("random-size,n", po::value<unsigned int>(), "randomly generated instance with n pairs")
        ("save,w", po::value<string>(), "save problem instance to (modifed) tsplib format file")
        ("solver-option,o", po::value<vector<string>>(), "solver option (e.g. foo=bar)")
        ("threads,p", po::value<unsigned int>(), "threads (default=1)")
        ("time-limit,t", po::value<unsigned int>(), "time limit in millis")
        ("soln-limit,l", po::value<unsigned int>(), "stop after n solutions are found")
        ;

    // Process command line options.
    po::variables_map varmap;
    try {
        po::store(po::command_line_parser(argc, argv).options(desc).run(), varmap);
    } catch (po::error& e) {
        cerr << e.what() << endl;
        cerr << desc << endl;
        return 1;
    }

    // Parse out solver options
    map<string, string> solver_options;
    if (varmap.count("solver-option") > 0)
        for (auto option_string : varmap["solver-option"].as<vector<string>>()) {
            // Option strings can't contain commas.
            if (option_string.find(",") != string::npos) {
                cerr << "invalid option string \"" << option_string << "\"" << endl;
                return 1;
            }

            // Split them on "=".
            vector<string> key_value;
            boost::algorithm::split(key_value, option_string, boost::is_any_of("="));
            if (key_value.size() != 2 || key_value[0] == "" || key_value[1] == "") {
                cerr << "invalid option string \"" << option_string << "\"" << endl;
                return 1;
            }

            solver_options[key_value[0]] = key_value[1];
        }

    if (varmap.count("help")) {
        cout << desc << endl;
        return 0;
    }

    // Choose a TSPPD solver.
    if (varmap.count("solver") != 1) {
        cerr << "solver required" << endl;
        return 1;
    }

    // Problem instance.
    TSPPD::Data::TSPPDProblem problem;

    if (varmap.count("input") > 0 && varmap.count("random-size") > 0) {
        cerr << "use either randomly generate input or a problem file" << endl;
        return 1;

    } else if (varmap.count("input") == 1) {
        problem = TSPPD::IO::TSPProblemReader::read(varmap["input"].as<string>());

    } else if (varmap.count("random-size") == 1) {
        // Size (number of pairs) must be > 0.
        unsigned size = varmap["random-size"].as<unsigned int>();
        if (size < 1) {
            cerr << "random-size must be > 0" << endl;
            return 1;
        }

        // Get a random seed.
        unsigned int seed = 0;
        if (varmap.count("random-seed") == 1)
            seed = varmap["random-seed"].as<unsigned int>();

        problem = TSPPD::Data::TSPPDProblemGenerator::generate(size, seed);

    } else {
        cerr << "input or random-size are required" << endl;
        return 1;
    }

    // Make sure the problem is valid.
    if (problem.nodes.size() < 2) {
        cerr << "invalid problem" << endl;
        return 1;
    }

    // Save the problem instance.
    if (varmap.count("save") == 1) {
        try {
            TSPPD::IO::TSPProblemWriter::write(varmap["save"].as<string>(), problem);
        } catch (TSPPD::Util::TSPPDException &e) {
            cerr << "error: " << e.what() << endl;
            return 1;
        }
    }

    // Name of the solver.
    auto solver_abbrev = varmap["solver"].as<string>();

    // Number of solver threads (MIP & CP only).
    unsigned int threads = 1;
    if (varmap.count("threads") == 1) {
        threads = varmap["threads"].as<unsigned int>();
        if (threads < 1)
            threads = 1;
    }

    // What is our output format?
    auto format = TSPPD::IO::HUMAN;
    if (varmap.count("format") == 1) {
        auto format_str = varmap["format"].as<string>();

        if (format_str == "csv")
            format = TSPPD::IO::CSV;
        else if (format_str == "human")
            format = TSPPD::IO::HUMAN;
        else {
            cerr << "invalid format: " << format_str << endl;
            return 1;
        }
    }

    TSPPD::IO::TSPSolutionWriter writer(problem, solver_abbrev, threads, solver_options, format);

    try {
        // Instantiate the solver.
        shared_ptr<TSPSolver> solver;

        if (solver_abbrev == "tsp-enum")
            solver = make_shared<EnumerativeTSPSolver>(problem, solver_options, writer);
        else if (solver_abbrev == "tsppd-enum")
            solver = make_shared<EnumerativeTSPPDSolver>(problem, solver_options, writer);

        else if (solver_abbrev == "tsp-focacci")
            solver = make_shared<FocacciTSPSolver>(problem, solver_options, writer);
        else if (solver_abbrev == "tsppd-focacci")
            solver = make_shared<FocacciTSPPDSolver>(problem, solver_options, writer);

        else if (solver_abbrev == "tsp-ruland")
            solver = make_shared<RulandTSPSolver>(problem, solver_options, writer);
        else if (solver_abbrev == "tsppd-ruland")
            solver = make_shared<RulandTSPPDSolver>(problem, solver_options, writer);
        else if (solver_abbrev == "tsppd-ruland+")
            solver = make_shared<RulandTSPPDPlusSolver>(problem, solver_options, writer);

        else if (solver_abbrev == "tsp-sarin")
            solver = make_shared<SarinTSPSolver>(problem, solver_options, writer);
        else if (solver_abbrev == "tsppd-sarin")
            solver = make_shared<SarinTSPPDSolver>(problem, solver_options, writer);

        else {
            cerr << "unknown solver: " << solver_abbrev << endl;
            return 1;
        }

        // Set a time limit if there is one
        if (varmap.count("time-limit") == 1)
            solver->time_limit = varmap["time-limit"].as<unsigned int>();

        // Set a solution limit if there is one
        if (varmap.count("soln-limit") == 1)
            solver->solution_limit = varmap["soln-limit"].as<unsigned int>();

        // Set a thread count
        solver->threads = threads;

        // Write a CSV output header.
        if (varmap.count("no-header") < 1)
            writer.write_header();

        solver->solve();

    } catch (GRBException &e) {
        cerr << "gurobi error: " << e.getMessage() << endl;
        return 1;

    } catch (TSPPD::Util::TSPPDException &e) {
        cerr << "error: " << e.what() << endl;
        return 1;
    }

    return 0;
}
