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

#include <fstream>

#include <tsppd/io/tsp_problem_writer.h>
#include <tsppd/util/exception.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Util;
using namespace std;

void TSPProblemWriter::write(const string filename, const TSPPDProblem& problem) {
    ofstream os(filename);

    if (os.is_open()) {
        os << "NAME: " << problem.name << endl
           << "TYPE: TSP" << endl
           << "COMMENT: " << problem.comment << endl
           << "DIMENSION: " << problem.nodes.size() << endl
           << "EDGE_WEIGHT_TYPE: EUC_2D" << endl
           << "NODE_COORD_SECTION" << endl;

        for (auto node : problem.nodes) {
            auto coord = problem.coordinate(node);
            os << node << " " << coord.first << " " << coord.second << endl;
        }

        auto pickups = problem.pickups();
        if (pickups.size() > 0) {
            os << "PRECEDENCE_SECTION" << endl;
            os << problem.nodes[0] << " " << problem.successor(problem.nodes[0]) << endl;
            for (auto p : pickups)
                os << p << " " << problem.successor(p) << endl;
        }

        os << "EOF" << endl;
        os.close();

    } else {
        throw TSPPDException("unable to write file '" + filename + "'");
    }
}
