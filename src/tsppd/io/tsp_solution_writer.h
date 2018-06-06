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

#ifndef TSPPD_IO_TSP_SOLUTION_WRITER_H
#define TSPPD_IO_TSP_SOLUTION_WRITER_H

#include <ctime>
#include <map>
#include <ostream>

#include <tsppd/data/tsppd_problem.h>
#include <tsppd/data/tsppd_search_statistics.h>
#include <tsppd/data/tsppd_solution.h>

namespace TSPPD {
    namespace IO {
        enum TSPSolutionFormat { HUMAN, CSV };
        const unsigned int TSPWriterSeparatorLength = 150;

        class TSPSolutionWriter {
        public:
            TSPSolutionWriter(
                const TSPPD::Data::TSPPDProblem& problem,
                const std::string solver,
                const unsigned int threads,
                const std::map<std::string, std::string> options,
                const TSPSolutionFormat format
            );

            void write_header();
            void write(const TSPPD::Data::TSPPDSearchStatistics& stats, const bool force = false);

       protected:
            const TSPPD::Data::TSPPDProblem problem;
            const std::string solver;
            const unsigned int threads;
            const std::map<std::string, std::string> options;
            const TSPSolutionFormat format;

            struct timespec start_wall;
            const clock_t start_cpu;

            // For avoiding duplicate output in human mode.
            std::string last_dual_str = "";
            std::string last_primal_str = "";
        };
    }
}

#endif
