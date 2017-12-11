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

#ifndef TSPPD_IO_TSP_PROBLEM_READER_H
#define TSPPD_IO_TSP_PROBLEM_READER_H

#include <fstream>
#include <vector>
#include <utility>

#include <tsppd/data/tsppd_problem.h>

namespace TSPPD {
    namespace IO {
        class TSPProblemReader {
        public:
            static TSPPD::Data::TSPPDProblem read(std::string filename);

        private:
            static std::vector<std::vector<int>> read_edge_weights(std::vector<std::string> lines);
            static std::pair<std::vector<std::string>, std::vector<std::pair<double, double>>>
                read_coordinate(std::vector<std::string> lines);
            static std::vector<std::pair<std::string, std::string>>
                read_precedence(std::vector<std::string> lines);

       };
    }
}

#endif
