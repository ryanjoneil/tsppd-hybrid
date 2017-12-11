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

#ifndef TSPPD_SOLVER_ENUMERATIVE_TSPPD_SOLVER_H
#define TSPPD_SOLVER_ENUMERATIVE_TSPPD_SOLVER_H

#include <tsppd/data/tsppd_arc.h>
#include <tsppd/solver/enumerative/enumerative_tsp_solver.h>

namespace TSPPD {
    namespace Solver {
        class EnumerativeTSPPDSolver : public EnumerativeTSPSolver {
        public:
            EnumerativeTSPPDSolver(
                const TSPPD::Data::TSPPDProblem& problem,
                const std::map<std::string, std::string> options,
                TSPPD::IO::TSPSolutionWriter& writer
            );

            virtual std::string name() const override { return "tsppd-enum"; }

        protected:
            bool feasible(TSPPD::Data::TSPPDArc next) override;
       };
    }
}

#endif
