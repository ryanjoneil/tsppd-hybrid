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

#ifndef TSPPD_SOLVER_GECODE_CLOSEST_NEIGHBOR_DUAL_H
#define TSPPD_SOLVER_GECODE_CLOSEST_NEIGHBOR_DUAL_H

#include <vector>

#include <gecode/int.hh>

#include <tsppd/data/tsppd_arc.h>
#include <tsppd/data/tsppd_problem.h>

namespace TSPPD {
    namespace Solver {
        class GecodeClosestNeighborDual : public Gecode::Propagator {
        public:
            GecodeClosestNeighborDual(
                Gecode::Space& home,
                Gecode::Int::IntView next,
                Gecode::Int::IntView closest_cost,
                const unsigned int node_index,
                const TSPPD::Data::TSPPDProblem& problem
            );

            GecodeClosestNeighborDual(
                Gecode::Space& home,
                bool share,
                GecodeClosestNeighborDual& p
            );

            virtual Gecode::Propagator* copy(Gecode::Space& home, bool share);
            virtual size_t dispose(Gecode::Space& home);

            virtual Gecode::PropCost cost(const Gecode::Space& home, const Gecode::ModEventDelta& med) const;
            virtual void reschedule(Gecode::Space& home);
            virtual Gecode::ExecStatus propagate(Gecode::Space& home, const Gecode::ModEventDelta& med);

            static Gecode::ExecStatus post(
                Gecode::Space& home,
                Gecode::Int::IntView next,
                Gecode::Int::IntView closest_cost,
                const unsigned int node_index,
                const TSPPD::Data::TSPPDProblem& problem
            );

        protected:
            Gecode::Int::IntView next;
            Gecode::Int::IntView closest_cost;
            const unsigned int node_index;
            const TSPPD::Data::TSPPDProblem& problem;
            unsigned int arc_index;
        };

        void closest_neighbor_dual(
            Gecode::Space& home,
            Gecode::IntVarArray next,
            Gecode::IntVar dual,
            const TSPPD::Data::TSPPDProblem& problem
        );
    }
}

#endif
