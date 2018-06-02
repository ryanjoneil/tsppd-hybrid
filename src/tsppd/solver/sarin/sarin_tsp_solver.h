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

#ifndef TSPPD_SOLVER_SARIN_TSP_SOLVER_H
#define TSPPD_SOLVER_SARIN_TSP_SOLVER_H

#include <map>
#include <memory>
#include <utility>

#include <gurobi_c++.h>

#include <tsppd/solver/tsp_solver.h>

namespace TSPPD {
    namespace Solver {
        // MIP TSPPD Solver based on:
        //
        // Subhash C. Sarin, Hanif D. Sherali, and Ajay Bhootra. "New tighter polynomial length formulations for
        // the asymmetric traveling salesman problem with and without precedence constraints." Operations research
        // letters 33, no. 1 (2005): 62-70.
        class SarinTSPSolver : public TSPSolver {
        public:
            SarinTSPSolver(
                const TSPPD::Data::TSPPDProblem& problem,
                const std::map<std::string, std::string> options,
                TSPPD::IO::TSPSolutionWriter& writer
            );

            std::string name() const { return "tsp-sarin"; }
            virtual TSPPD::Data::TSPPDSolution solve();

        protected:
            void initialize_variables();
            void initialize_constraints();
            std::vector<unsigned int> get_path();

            GRBEnv env;
            GRBModel model;
            std::vector<std::vector<GRBVar>> x;
            std::vector<std::vector<GRBVar>> y;
       };
    }
}

#endif
