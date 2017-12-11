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

#ifndef TSPPD_SOLVER_GECODE_TSP_BRANCHER_H
#define TSPPD_SOLVER_GECODE_TSP_BRANCHER_H

#include <map>
#include <ostream>
#include <vector>

#include <gecode/int.hh>

#include <tsppd/data/tsppd_arc.h>
#include <tsppd/data/tsppd_problem.h>
#include <tsppd/solver/gecode/brancher/gecode_tsp_branch_choice.h>

namespace TSPPD {
    namespace Solver {
        enum GecodeTSPBrancherType { BRANCHER_CN, BRANCHER_SEQ_CN, BRANCHER_REGRET };

        class GecodeTSPBrancher : public Gecode::Brancher {
        public:
            GecodeTSPBrancher(
                Gecode::Home home,
                Gecode::ViewArray<Gecode::Int::IntView>& next,
                const TSPPD::Data::TSPPDProblem& problem
            );

            GecodeTSPBrancher(Gecode::Space& home, bool share, GecodeTSPBrancher& b);
            virtual Gecode::Actor* copy(Gecode::Space& home, bool share) = 0;
            virtual size_t dispose(Gecode::Space& home);

            virtual bool status(const Gecode::Space& home) const;
            virtual Gecode::Choice* choice(Gecode::Space& home) = 0;
            virtual Gecode::Choice* choice(const Gecode::Space&, Gecode::Archive& e);
            virtual Gecode::ExecStatus commit(Gecode::Space& home, const Gecode::Choice& c, unsigned int a);

            virtual void print(
                const Gecode::Space& home,
                const Gecode::Choice& c,
                unsigned int a,
                std::ostream& o
            ) const;

        protected:
            unsigned int closest_feasible_arc_index(unsigned int from, unsigned int start);

            Gecode::ViewArray<Gecode::Int::IntView> next;
            const TSPPD::Data::TSPPDProblem& problem;
            std::vector<unsigned int> indexes;
        };
    }
}

#endif
