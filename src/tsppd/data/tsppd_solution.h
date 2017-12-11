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

#ifndef TSPPD_DATA_TSPPD_SOLUTION_H
#define TSPPD_DATA_TSPPD_SOLUTION_H

#include <vector>

#include <tsppd/data/tsppd_problem.h>

namespace TSPPD {
    namespace Data {
        class TSPPDSolution {
        public:
            TSPPDSolution(const TSPPD::Data::TSPPDProblem& problem, const std::vector<std::string> order);
            TSPPDSolution(const TSPPD::Data::TSPPDProblem& problem, const std::vector<unsigned int> index_order);

            bool feasible() const;

            const TSPPD::Data::TSPPDProblem& problem;
            const std::vector<std::string> order;
            const int cost;

       private:
            std::vector<std::string> string_order(std::vector<unsigned int> index_order);
            int compute_cost();
       };
    }
}

#endif
