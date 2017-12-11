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
#include <iostream>
#include <memory>
#include <vector>

#include <tsppd/ap/ap_solver.h>
#include <tsppd/ap/gurobi_ap_solver.h>
#include <tsppd/ap/primal_dual_ap_solver.h>

using namespace TSPPD::AP;
using namespace std;

int main(int argc, char** argv) {
    vector<vector<int>> obj = {
        {25, 40,  9, 15},
        {12, 33, 22, 31},
        {18, 37, 18, 30},
        {32, 11, 12, 26}
    };

    cout << "AP Solver: Gurobi" << endl;
    GurobiAPSolver s(obj.size());

    for (unsigned int row = 0; row < obj.size(); ++row)
        for (unsigned int col = 0; col < obj[row].size(); ++col)
            s.set_obj({row, col}, obj[row][col]);

    s.solve();
    s.print(cout);
    s.set_bounds({2, 2}, false, false);
    s.set_bounds({1, 0}, false, false);
    s.set_bounds({0, 3}, false, false);
    s.set_bounds({3, 1}, false, false);
    s.set_bounds({2, 0}, false, false);
    s.set_bounds({3, 3}, false, false);
    s.solve();
    s.print(cout);
    cout << endl;

    cout << "AP Solver: Primal-Dual" << endl;

    PrimalDualAPSolver h(obj.size());

    for (unsigned int row = 0; row < obj.size(); ++row)
        for (unsigned int col = 0; col < obj[row].size(); ++col)
            h.set_obj({row, col}, obj[row][col]);

    h.solve();
    h.print(cout);
    h.set_bounds({2, 2}, false, false);
    h.set_bounds({1, 0}, false, false);
    h.set_bounds({0, 3}, false, false);
    h.set_bounds({3, 1}, false, false);
    h.set_bounds({2, 0}, false, false);
    h.set_bounds({3, 3}, false, false);
    h.solve();
    h.print(cout);
}
