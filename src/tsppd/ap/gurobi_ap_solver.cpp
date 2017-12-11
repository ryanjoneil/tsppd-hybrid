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

#include <cmath>

#include <tsppd/ap/gurobi_ap_solver.h>

using namespace TSPPD::AP;
using namespace std;

GurobiAPSolver::GurobiAPSolver(const unsigned int size) :
    APSolver(size),
    env(),
    model(env) {

    initialize_model();
}

bool GurobiAPSolver::solve() {
    model.optimize();

    if (model.get(GRB_IntAttr_Status) != GRB_OPTIMAL)
        return false;

    return true;
}

void GurobiAPSolver::set_obj(APIndex index, int obj) {
    x[index.first][index.second].set(GRB_DoubleAttr_Obj, obj);
}

void GurobiAPSolver::set_bounds(APIndex index, bool lb, bool ub) {
    x[index.first][index.second].set(GRB_DoubleAttr_LB, lb);
    x[index.first][index.second].set(GRB_DoubleAttr_UB, ub);
}

int GurobiAPSolver::get_z() const {
    return model.get(GRB_DoubleAttr_ObjVal);
}

bool GurobiAPSolver::get_x(APIndex index) {
    return *model.get(GRB_DoubleAttr_X, &(x[index.first][index.second]), 1) > 0.5;
}

int GurobiAPSolver::get_rc(APIndex index) {
    return round(*model.get(GRB_DoubleAttr_RC, &(x[index.first][index.second]), 1));
}

bool GurobiAPSolver::get_lb(APIndex index) {
    return *model.get(GRB_DoubleAttr_LB, &(x[index.first][index.second]), 1) > 0.5;
}

bool GurobiAPSolver::get_ub(APIndex index) {
    return *model.get(GRB_DoubleAttr_UB, &(x[index.first][index.second]), 1) > 0.5;
}

int GurobiAPSolver::get_u(unsigned int row) {
    return round(*model.get(GRB_DoubleAttr_Pi, &(u[row]), 1));
}

int GurobiAPSolver::get_v(unsigned int col) {
    return round(*model.get(GRB_DoubleAttr_Pi, &(v[col]), 1));
}

void GurobiAPSolver::initialize_model() {
    // Silence output.
    model.getEnv().set(GRB_IntParam_OutputFlag, 0);

    // x = AP vars
    for (unsigned int row = 0; row < size; ++row) {
        vector<GRBVar> x_i;

        for (unsigned int col = 0; col < size; ++col)
            x_i.push_back(model.addVar(0, 1, 0, GRB_CONTINUOUS));

        x.push_back(x_i);
    }

    // Constraints
    for (unsigned int row = 0; row < size; ++row) {
        GRBLinExpr expr_row = 0;
        GRBLinExpr expr_col = 0;

        for (unsigned int col = 0; col < size; ++col) {
            expr_row += x[row][col];
            expr_col += x[col][row];
        }

        u.push_back(model.addConstr(expr_row == 1));
        v.push_back(model.addConstr(expr_col == 1));
    }
}
