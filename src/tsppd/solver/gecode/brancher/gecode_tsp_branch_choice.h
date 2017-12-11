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

#ifndef TSPPD_SOLVER_GECODE_BRANCH_CHOICE_H
#define TSPPD_SOLVER_GECODE_BRANCH_CHOICE_H

#include <gecode/int.hh>

#include <tsppd/data/tsppd_problem.h>

namespace TSPPD {
    namespace Solver {
        class GecodeTSPBranchChoice : public Gecode::Choice {
        public:
            GecodeTSPBranchChoice(const Gecode::Brancher& b, int index, int value) :
                Gecode::Choice(b, 2), next_index(index), next_value(value) { }

            virtual size_t size(void) const {
                return sizeof(*this);
            }

            virtual void archive(Gecode::Archive& e) const {
                Gecode::Choice::archive(e);
                e << next_index << next_value;
            }

            const int next_index;
            const int next_value;
        };
    }
}

#endif
