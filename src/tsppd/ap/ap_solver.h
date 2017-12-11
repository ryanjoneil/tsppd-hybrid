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

#ifndef TSPPD_AP_SOLVER_H
#define TSPPD_AP_SOLVER_H

#include <ostream>

namespace TSPPD {
    namespace AP {
        typedef std::pair<unsigned int, unsigned int> APIndex;

        class APSolver {
        public:
            APSolver(const unsigned int size);
            virtual ~APSolver() { }

            virtual bool solve() = 0;

            virtual void set_obj(APIndex index, int obj) = 0;
            virtual void set_bounds(APIndex index, bool lb, bool ub) = 0;

            virtual int get_z() const = 0;
            virtual bool get_x(APIndex index) = 0;
            virtual bool get_lb(APIndex index) = 0;
            virtual bool get_ub(APIndex index) = 0;
            virtual int get_rc(APIndex index) = 0;
            virtual int get_u(unsigned int row) = 0;
            virtual int get_v(unsigned int col) = 0;

            void print(std::ostream& out);

        protected:
            const unsigned int size;
        };
    }
}

#endif
