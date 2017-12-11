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

#ifndef TSPPD_PRIMAL_DUAL_AP_SOLVER_H
#define TSPPD_PRIMAL_DUAL_AP_SOLVER_H

#include <set>
#include <vector>

#include <tsppd/ap/ap_solver.h>

// This is a C++ port of:
//
// Carpaneto, Giorgio, Silvano Martello, and Paolo Toth.
// "Algorithms and codes for the assignment problem."
// Annals of operations research 13, no. 1 (1988): 191-223.

namespace TSPPD {
    namespace AP {
        class PrimalDualAPSolver : public APSolver {
        public:
            PrimalDualAPSolver(const unsigned int size);

            virtual bool solve() override;

            virtual void set_obj(APIndex index, int _obj) override;
            virtual void set_bounds(APIndex index, bool _lb, bool _ub) override;

            virtual int get_z() const override;
            virtual bool get_x(APIndex index) override;
            virtual bool get_lb(APIndex index) override;
            virtual bool get_ub(APIndex index) override;
            virtual int get_rc(APIndex index) override;
            virtual int get_u(unsigned int row) override;
            virtual int get_v(unsigned int col) override;

        protected:
            int z;
            std::vector<std::vector<int>> a;
            std::vector<std::vector<bool>> lb;
            std::vector<std::vector<bool>> ub;
            std::vector<int> u;
            std::vector<int> v;

            std::vector<APIndex> zeroes;
            std::vector<bool> covered_rows;
            std::vector<bool> covered_cols;

        private:
            void initialize();
            int path(int i);
            void increase(int i, int j);

            void initialize_phase_1();
            void initialize_phase_2();

            int min_row(int j);
            int min_col(int i);

            std::vector<int> f;        // f[i] = j if i assigned to j, -1 if i unassigned
            std::vector<int> f_bar;    // f_bar[j] = row assigned to column j, -1 if unassigned
            std::vector<int> p;        // first column of row i not yet examined in phase 2

            std::vector<int> LR;       // vector of labelled rows
            std::set<int> UC;          // set of unlabelled columns
            std::vector<int> c;        // c[j] = row preceding column j in current alternating path
            std::vector<int> pi;       // pi[j] = min { a[i,j] - u[i] - v[j] | i in LR, i != f_bar[j] }

            bool initialized;
        };
    }
}

#endif
