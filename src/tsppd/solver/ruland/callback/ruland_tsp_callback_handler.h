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

#ifndef TSPPD_SOLVER_RULAND_TSP_CALLBACK_HANDLER_H
#define TSPPD_SOLVER_RULAND_TSP_CALLBACK_HANDLER_H

#include <map>
#include <utility>
#include <vector>

#include <gurobi_c++.h>

#include <tsppd/data/tsppd_problem.h>
#include <tsppd/io/tsp_solution_writer.h>
#include <tsppd/solver/ruland/callback/ruland_tsp_callback.h>
#include <tsppd/solver/ruland/ruland_subtour_finder.h>
#include <tsppd/solver/tsp_solver.h>

namespace TSPPD {
    namespace Solver {
        class RulandTSPCallbackHandler : public GRBCallback {
        public:
            RulandTSPCallbackHandler(
                const TSPPD::Solver::TSPSolver* solver,
                const TSPPD::Data::TSPPDProblem& problem,
                std::map<std::pair<unsigned int, unsigned int>, GRBVar> arcs,
                std::vector<std::shared_ptr<RulandTSPCallback>> callbacks,
                TSPPD::IO::TSPSolutionWriter& writer
            );

            // These wrap callback methods so we can break our callbacks into multiple classes.
            double get_solution(GRBVar v);
            double* get_solution(const GRBVar* xvars, int len);
            void add_lazy(const GRBTempConstr& tc);
            void add_lazy(const GRBLinExpr& expr, char sense, double rhs);

        protected:
            void callback();

            const TSPPD::Solver::TSPSolver* solver;
            const TSPPD::Data::TSPPDProblem& problem;
            std::map<std::pair<unsigned int, unsigned int>, GRBVar> arcs;
            std::vector<std::shared_ptr<RulandTSPCallback>> callbacks;
            TSPPD::IO::TSPSolutionWriter writer;
            TSPPD::Solver::RulandSubtourFinder subtour_finder;
       };
    }
}

#endif
