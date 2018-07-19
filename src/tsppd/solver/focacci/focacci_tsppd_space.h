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

#ifndef TSPPD_SOLVER_FOCACCI_TSPPD_SPACE_H
#define TSPPD_SOLVER_FOCACCI_TSPPD_SPACE_H

#include <vector>

#include <gecode/set.hh>
#include <gecode/minimodel.hh>

#include <tsppd/data/tsppd_problem.h>
#include <tsppd/solver/focacci/propagator/focacci_tsppd_precede_propagator.h>
#include <tsppd/solver/focacci/focacci_tsp_space.h>

namespace TSPPD {
    namespace Solver {
        class FocacciTSPPDSpace : public FocacciTSPSpace {
        public:
            FocacciTSPPDSpace(const TSPPD::Data::TSPPDProblem& problem);

            // Search & optimization support
            FocacciTSPPDSpace(FocacciTSPPDSpace& rs);
            virtual Gecode::Space* copy();

            virtual void initialize_constraints();
            void initialize_precedence_propagators(const FocacciTSPPDPrecedePropagatorType precede_type);
            void initialize_additive_propagator();
            void initialize_assignment_propagator();
            void initialize_heldkarp_propagator();
            void initialize_omc_constraints();
        };
    }
}

#endif
