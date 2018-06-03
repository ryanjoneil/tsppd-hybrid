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

#ifndef TSPPD_SOLVER_RULAND_SUBTOUR_FINDER_H
#define TSPPD_SOLVER_RULAND_SUBTOUR_FINDER_H

#include <map>
#include <set>
#include <utility>
#include <vector>

#include <tsppd/data/tsppd_problem.h>
#include <tsppd/io/tsp_solution_writer.h>
#include <tsppd/solver/ruland/ruland_subtour_finder.h>

namespace TSPPD {
    namespace Solver {
        class RulandSubtourFinder {
        public:
            RulandSubtourFinder(const TSPPD::Data::TSPPDProblem& problem);

            std::vector<std::vector<unsigned int>> subtours(
                const std::map<std::pair<unsigned int, unsigned int>, bool>& arcs
            );

        protected:
            const TSPPD::Data::TSPPDProblem& problem;

        private:
            std::map<unsigned int, std::set<unsigned int>> connections(
                const std::map<std::pair<unsigned int, unsigned int>, bool>& arcs
            );
       };
    }
}

#endif
