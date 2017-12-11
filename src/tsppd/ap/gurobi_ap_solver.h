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

#ifndef TSPPD_GUROBI_AP_SOLVER_H
#define TSPPD_GUROBI_AP_SOLVER_H

#include <vector>

#include <gurobi_c++.h>

#include <tsppd/ap/ap_solver.h>

namespace TSPPD {
    namespace AP {
        class GurobiAPSolver : public APSolver {
        public:
            GurobiAPSolver(const unsigned int size);

            virtual bool solve() override;

            virtual void set_obj(APIndex index, int obj) override;
            virtual void set_bounds(APIndex index, bool lb, bool ub) override;

            virtual int get_z() const override;
            virtual bool get_x(APIndex index) override;
            virtual bool get_lb(APIndex index) override;
            virtual bool get_ub(APIndex index) override;
            virtual int get_rc(APIndex index) override;
            virtual int get_u(unsigned int row) override;
            virtual int get_v(unsigned int col) override;

        protected:
            void initialize_model();

            GRBEnv env;
            GRBModel model;
            std::vector<std::vector<GRBVar>> x;
            std::vector<GRBConstr> u;
            std::vector<GRBConstr> v;

            bool solved = false;
        };
    }
}

#endif
