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

#ifndef TSPPD_SOLVER_FOCACCI_TSP_HKAP_FILTER_H
#define TSPPD_SOLVER_FOCACCI_TSP_HKAP_FILTER_H

#include <vector>

#include <gecode/int.hh>
#include <gecode/minimodel.hh>

#include <tsppd/data/tsppd_problem.h>
#include <tsppd/solver/focacci/filter/focacci_tsp_heldkarp_filter.h>

namespace TSPPD {
    namespace Solver {
        class FocacciTSPHKAPFilter : public FocacciTSPHeldKarpFilter {
        public:
            FocacciTSPHKAPFilter(
                Gecode::Home home,
                Gecode::ViewArray<Gecode::Int::IntView>& next,
                Gecode::Int::IntView& primal,
                const TSPPD::Data::TSPPDProblem& problem,
                const unsigned int max_iterations
            );

            FocacciTSPHKAPFilter(Gecode::Space& home, FocacciTSPHKAPFilter& p);

            virtual Gecode::Propagator* copy(Gecode::Space& home);
            virtual Gecode::ExecStatus propagate(Gecode::Space& home, const Gecode::ModEventDelta& med);

            static Gecode::ExecStatus post(
                Gecode::Home home,
                Gecode::ViewArray<Gecode::Int::IntView>& next,
                Gecode::Int::IntView& primal,
                const TSPPD::Data::TSPPDProblem& problem,
                const unsigned int max_iterations
            );
        };

        void tsppd_hkap(
            Gecode::Home home,
            Gecode::IntVarArray& next,
            Gecode::IntVar& primal,
            const TSPPD::Data::TSPPDProblem& problem,
            const unsigned int max_iterations
        );
    }
}

#endif
