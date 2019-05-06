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

#include <iomanip>

#include <tsppd/ap/ap_solver.h>

using namespace TSPPD::AP;
using namespace std;

APSolver::APSolver(const unsigned int size) : size(size) { }

void APSolver::print(ostream& out) {
    out << "========================================================================" << endl;
    out << "z\t" << setw(5) << setfill(' ') << right << get_z() << endl << endl;

    out << "u\t";
    for (unsigned int row = 0; row < size; ++row)
        out << setw(5) << setfill(' ') << right << get_u(row);
    out << endl;

    out << "v\t";
    for (unsigned int col = 0; col < size; ++col)
        out << setw(5) << setfill(' ') << right << get_v(col);
    out << endl << endl;

    out << "bnd\t";
    for (unsigned int row = 0; row < size; ++row) {
        for (unsigned int col = 0; col < size; ++col) {
            auto lb = get_lb({row, col});
            auto ub = get_ub({row, col});

            out << setw(5) << setfill(' ') << right;

            if (!ub)
                out << "0";
            else if (lb)
                out << "1";
            else
                out << "-";
        }
        out << endl << "\t";
    }
    out << endl;

    out << "rc\t";
    for (unsigned int row = 0; row < size; ++row) {
        for (unsigned int col = 0; col < size; ++col) {
            out << setw(5) << setfill(' ') << right;
            if (get_ub({row, col}))
                out << get_rc({row, col});
            else
                out << "-";
        }
        out << endl << "\t";
    }
    out << endl;

    out << "x\t";
    for (unsigned int row = 0; row < size; ++row) {
        for (unsigned int col = 0; col < size; ++col) {
            out << setw(5) << setfill(' ') << right
                << (get_x({row, col}) ? "1" : "-");
        }
        out << endl;
        if (row < size - 1)
            out << "\t";
    }

    out << "------------------------------------------------------------------------" << endl;
}
