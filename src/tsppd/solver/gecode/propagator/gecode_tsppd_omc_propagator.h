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

#ifndef TSPPD_SOLVER_GECODE_TSPPD_PROPAGATOR_H
#define TSPPD_SOLVER_GECODE_TSPPD_PROPAGATOR_H

#include <gecode/int.hh>

#include <tsppd/data/tsppd_problem.h>

namespace TSPPD {
    namespace Solver {
        class GecodeTSPPDOMCPropagator : public Gecode::Propagator {
        public:
            GecodeTSPPDOMCPropagator(
                Gecode::Space& home,
                Gecode::ViewArray<Gecode::Int::IntView>& next,
                const unsigned int index,
                const TSPPD::Data::TSPPDProblem& problem
            );

            GecodeTSPPDOMCPropagator(Gecode::Space& home, GecodeTSPPDOMCPropagator& p);

            virtual Gecode::Propagator* copy(Gecode::Space& home);
            virtual size_t dispose(Gecode::Space& home);

            virtual Gecode::PropCost cost(const Gecode::Space& home, const Gecode::ModEventDelta& med) const;
            virtual void reschedule(Gecode::Space& home);
            virtual Gecode::ExecStatus propagate(Gecode::Space& home, const Gecode::ModEventDelta& med);

            static Gecode::ExecStatus post(
                Gecode::Space& home,
                Gecode::ViewArray<Gecode::Int::IntView>& next,
                const unsigned int index,
                const TSPPD::Data::TSPPDProblem& problem
            );

        protected:
            Gecode::ViewArray<Gecode::Int::IntView> next;
            const unsigned int index;
            const TSPPD::Data::TSPPDProblem& problem;

            const unsigned int start_index;
            const unsigned int end_index;
        };

        void tsppd_omc(
            Gecode::Space& home,
            Gecode::IntVarArray& next,
            const TSPPD::Data::TSPPDProblem& problem
        );
    }
}

#endif