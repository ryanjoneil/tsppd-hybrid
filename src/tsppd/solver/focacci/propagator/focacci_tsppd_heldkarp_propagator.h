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

#ifndef TSPPD_SOLVER_FOCACCI_TSPPD_HELDKARP_PROPAGATOR_H
#define TSPPD_SOLVER_FOCACCI_TSPPD_HELDKARP_PROPAGATOR_H

#include <set>
#include <utility>
#include <vector>

#include <gecode/int.hh>
#include <gecode/minimodel.hh>

#include <tsppd/ap/primal_dual_ap_solver.h>
#include <tsppd/data/tsppd_problem.h>

namespace TSPPD {
    namespace Solver {
        class FocacciTSPPDHeldKarpPropagator : public Gecode::Propagator {
        public:
            FocacciTSPPDHeldKarpPropagator(
                Gecode::Home home,
                Gecode::ViewArray<Gecode::Int::IntView>& next,
                Gecode::Int::IntView& primal,
                const TSPPD::Data::TSPPDProblem& problem
            );

            FocacciTSPPDHeldKarpPropagator(Gecode::Space& home, FocacciTSPPDHeldKarpPropagator& p);

            virtual Gecode::Propagator* copy(Gecode::Space& home);
            virtual size_t dispose(Gecode::Space& home);

            virtual Gecode::PropCost cost(const Gecode::Space& home, const Gecode::ModEventDelta& med) const;
            virtual void reschedule(Gecode::Space& home);
            virtual Gecode::ExecStatus propagate(Gecode::Space& home, const Gecode::ModEventDelta& med);

            static Gecode::ExecStatus post(
                Gecode::Home home,
                Gecode::ViewArray<Gecode::Int::IntView>& next,
                Gecode::Int::IntView& primal,
                const TSPPD::Data::TSPPDProblem& problem
            );

        protected:
            double one_tree(const std::vector<double>& potentials, std::vector<std::set<int>>& edges);
            int undirected_cost(int i, int j);
            double transformed_cost(int i, int j, const std::vector<double> potentials);
            int marginal_cost(
                int from,
                int to,
                const std::vector<std::set<int>>& edges
            );
            int marginal_cost(
                int from,
                int to,
                const std::vector<std::set<int>>& edges,
                std::vector<bool> seen,
                int node,
                int max_edge_cost
            );

            Gecode::ViewArray<Gecode::Int::IntView> next;
            Gecode::Int::IntView primal;
            const TSPPD::Data::TSPPDProblem& problem;

            const int start_index;
            const int end_index;

            const double C = 1.5;
        };

        void tsppd_heldkarp(
            Gecode::Home home,
            Gecode::IntVarArray& next,
            Gecode::IntVar& primal,
            const TSPPD::Data::TSPPDProblem& problem
        );
    }
}

#endif
