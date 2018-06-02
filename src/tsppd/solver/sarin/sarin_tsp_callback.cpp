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

#include <tsppd/data/tsppd_search_statistics.h>
#include <tsppd/data/tsppd_solution.h>
#include <tsppd/solver/sarin/sarin_tsp_callback.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Solver;
using namespace std;

SarinTSPCallback::SarinTSPCallback(
    const TSPPDProblem& problem,
    vector<vector<GRBVar>> x,
    TSPSolutionWriter& writer) :
    problem(problem), x(x), writer(writer) { }

void SarinTSPCallback::callback() {
    if (where == GRB_CB_MIP) {
        // Just log bounds.
        TSPPDSearchStatistics stats;
        stats.primal = getDoubleInfo(GRB_CB_MIP_OBJBST);
        stats.dual = max(0.0, getDoubleInfo(GRB_CB_MIP_OBJBND));
        writer.write(stats);
   
    } else if (where == GRB_CB_MIPSOL) {
        // Pull out a new solution.
        vector<unsigned int> path{0};
        unsigned int from = 0;
        for (unsigned int i = 0; i < problem.nodes.size() - 1; ++i) {
            for (unsigned int to = 0; to < problem.nodes.size(); ++to) {
                if (getSolution(x[from][to]) > 0.5) {
                    path.push_back(to);
                    from = to;
                    break;
                }
            }
        }

        TSPPDSolution solution(problem, path);
        TSPPDSearchStatistics stats(solution);
        stats.dual = max(0.0, getDoubleInfo(GRB_CB_MIPSOL_OBJBND));
        writer.write(stats);
    }
}
