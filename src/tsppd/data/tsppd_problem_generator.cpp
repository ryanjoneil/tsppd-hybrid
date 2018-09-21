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

#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

#include <tsppd/data/tsppd_problem_generator.h>

#include <iostream>

using namespace TSPPD::Data;
using namespace std;

TSPPDProblem TSPPDProblemGenerator::generate(unsigned int size, unsigned int seed) {
    string name = string("random-") + to_string(size) + "-" + to_string(seed);
    string comment = string("size=") + to_string(size) + " seed=" + to_string(seed);

    vector<string> nodes;
    vector<pair<double, double>> coordinates;
    vector<pair<string, string>> pickup_delivery_pairs;
    for (unsigned int i = 0; i <= size; ++i) {
        string pickup = string("+") + to_string(i);
        string delivery = string("-") + to_string(i);

        nodes.push_back(pickup);
        nodes.push_back(delivery);

        auto xp = rand() % 1000;
        auto yp = rand() % 1000;

        auto xd = rand() % 1000;
        auto yd = rand() % 1000;
        if (i == 0) {
             xd = xp;
             yd = yp;
        }

        coordinates.push_back({xp, yp});
        coordinates.push_back({xd, yd});

        pickup_delivery_pairs.push_back({pickup, delivery});
    }

    return TSPPDProblem(name, TSP, comment, size, EUC_2D, nodes, coordinates, pickup_delivery_pairs);
}
