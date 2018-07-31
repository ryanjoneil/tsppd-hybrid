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

#ifndef TSPPD_SOLVER_FOCACCI_TSP_FILTER_ONE_TREE_H
#define TSPPD_SOLVER_FOCACCI_TSP_FILTER_ONE_TREE_H

#include <set>
#include <utility>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include <gecode/int.hh>
#include <gecode/minimodel.hh>

#include <tsppd/ap/primal_dual_ap_solver.h>
#include <tsppd/data/tsppd_problem.h>

#define FOCACCI_TSP_FILTER_ONE_TREE_MAX_ITERATIONS 100

namespace TSPPD {
    namespace Solver {
        typedef boost::adjacency_list<
            boost::vecS,
            boost::vecS,
            boost::directedS,
            boost::no_property,
            boost::property<boost::edge_weight_t, int>
        > OneTreeGraph;
        typedef boost::graph_traits <OneTreeGraph>::edge_descriptor OneTreeEdge;
        typedef boost::property_map<OneTreeGraph, boost::edge_weight_t>::type OneTreeWeights;

        class OneTree {
        public:
            OneTree(
                Gecode::ViewArray<Gecode::Int::IntView>& next,
                const TSPPD::Data::TSPPDProblem& problem,
                const unsigned int max_iterations
            );

            OneTree(
                Gecode::ViewArray<Gecode::Int::IntView>& next,
                const TSPPD::Data::TSPPDProblem& problem,
                const unsigned int max_iterations,
                TSPPD::AP::PrimalDualAPSolver* ap
            );

            bool is_done();
            double bound();
            double improve();
            bool has_edge(int from, int to);
            int marginal_cost(int from, int to);

        protected:
            void initialize_one_tree();
            double minimize_one_tree();
            void update_one_tree();

            int undirected_cost(int i, int j);
            double transformed_cost(int i, int j);

            int marginal_cost(
                int from,
                int to,
                std::vector<bool> seen,
                int node,
                int max_edge_cost
            );

            int cost(int i, int j);

            Gecode::ViewArray<Gecode::Int::IntView> next;
            const TSPPD::Data::TSPPDProblem& problem;
            const unsigned int max_iterations;
            TSPPD::AP::PrimalDualAPSolver* ap;

            OneTreeGraph graph;
            OneTreeWeights weights;
            std::vector<double> potentials;
            std::vector<std::set<int>> edges;

            const int start_index;
            const int end_index;

            unsigned int iteration;
            double w;
            double t1, ti;
            bool done = false;

            const double EPSILON = 10e-7;
        };
    }
}

#endif
