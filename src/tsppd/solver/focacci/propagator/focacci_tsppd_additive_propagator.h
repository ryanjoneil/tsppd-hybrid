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

#ifndef TSPPD_SOLVER_FOCACCI_TSPPD_ADDITIVE_PROPAGATOR_H
#define TSPPD_SOLVER_FOCACCI_TSPPD_ADDITIVE_PROPAGATOR_H

#include <vector>

#include <gecode/int.hh>
#include <gecode/minimodel.hh>

#include <tsppd/data/tsppd_problem.h>
#include <tsppd/solver/focacci/one_tree/focacci_tsppd_one_tree.h>
#include <tsppd/solver/focacci/propagator/focacci_tsppd_assignment_propagator.h>

namespace TSPPD {
    namespace Solver {
        class FocacciTSPPDAdditivePropagator : public FocacciTSPPDAssignmentPropagator {
        public:
            FocacciTSPPDAdditivePropagator(
                Gecode::Home home,
                Gecode::ViewArray<Gecode::Int::IntView>& next,
                Gecode::Int::IntView& primal,
                const TSPPD::Data::TSPPDProblem& problem
            );

            FocacciTSPPDAdditivePropagator(Gecode::Space& home, FocacciTSPPDAdditivePropagator& p);

            virtual Gecode::Propagator* copy(Gecode::Space& home);
            virtual Gecode::ExecStatus propagate(Gecode::Space& home, const Gecode::ModEventDelta& med);

            static Gecode::ExecStatus post(
                Gecode::Home home,
                Gecode::ViewArray<Gecode::Int::IntView>& next,
                Gecode::Int::IntView& primal,
                const TSPPD::Data::TSPPDProblem& problem
            );
        };

        void tsppd_additive(
            Gecode::Home home,
            Gecode::IntVarArray& next,
            Gecode::IntVar& primal,
            const TSPPD::Data::TSPPDProblem& problem
        );
    }
}

#endif
