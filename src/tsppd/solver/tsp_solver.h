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

#ifndef TSPPD_SOLVER_TSP_SOLVER_H
#define TSPPD_SOLVER_TSP_SOLVER_H

#include <ctime>
#include <iostream>
#include <map>
#include <string>

#include <tsppd/data/tsppd_problem.h>
#include <tsppd/data/tsppd_solution.h>
#include <tsppd/io/tsp_solution_writer.h>

namespace TSPPD {
    namespace Solver {
        class TSPSolver {
        public:
            TSPSolver(
                const TSPPD::Data::TSPPDProblem& problem,
                const std::map<std::string, std::string> options,
                TSPPD::IO::TSPSolutionWriter& writer
            ) : problem(problem), options(options), writer(writer) { }

            virtual std::string name() const = 0;
            virtual TSPPD::Data::TSPPDSolution solve() = 0;

            unsigned int time_limit = 0;
            unsigned int solution_limit = 0;

        protected:
            const TSPPD::Data::TSPPDProblem& problem;
            std::map<std::string, std::string> options;
            TSPPD::IO::TSPSolutionWriter& writer;
       };
    }
}

#endif
