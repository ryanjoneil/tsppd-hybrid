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

#ifndef TSPPD_SOLVER_ONEIL_ATSPPD_SOLVER_H
#define TSPPD_SOLVER_ONEIL_ATSPPD_SOLVER_H

#include <map>
#include <memory>
#include <utility>

#include <gurobi_c++.h>

#include <tsppd/solver/ap/ap_atsppd_solver.h>

namespace TSPPD {
    namespace Solver {
        // A new MIP TSPPD solver model
        //
        // Solver Options:
        //     relax:  relax model and add SEC as violated {on|off} (default=off)
        //     sec:    relaxed SEC form that uses x variables {cutset|subtour} (default=subtour)
        class ONeilATSPPDSolver : public APATSPPDSolver {
        public:
            ONeilATSPPDSolver(
                const TSPPD::Data::TSPPDProblem& problem,
                const std::map<std::string, std::string> options,
                TSPPD::IO::TSPSolutionWriter& writer
            );

            std::string name() const { return "atsppd-oneil"; }
            virtual TSPPD::Data::TSPPDSolution solve();

        protected:
            void initialize_oneil_options();
            void initialize_oneil_variables();
            void initialize_oneil_constraints();

            GRBLinExpr s(unsigned int i, unsigned int j);

            std::map<std::pair<unsigned int, unsigned int>, std::vector<GRBVar>> w;
       };
    }
}

#endif
