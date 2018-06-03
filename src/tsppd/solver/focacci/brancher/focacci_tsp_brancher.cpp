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

#include <tsppd/solver/focacci/brancher/focacci_tsp_brancher.h>

using namespace Gecode;
using namespace TSPPD::Data;
using namespace TSPPD::Solver;
using namespace std;

FocacciTSPBrancher::FocacciTSPBrancher(
    Home home,
    ViewArray<Int::IntView>& next,
    const TSPPDProblem& problem) :
    Brancher(home),
    next(next),
    problem(problem),
    indexes(vector<unsigned int>(next.size(), 0)) { }

FocacciTSPBrancher::FocacciTSPBrancher(Space& home, FocacciTSPBrancher& b) :
    Brancher(home, b),
    next(b.next),
    problem(b.problem),
    indexes(b.indexes) {

    next.update(home, b.next);
}

size_t FocacciTSPBrancher::dispose(Gecode::Space& home) {
    indexes.~vector<unsigned int>();
    (void) Brancher::dispose(home);
    return sizeof(*this);
}

bool FocacciTSPBrancher::status(const Space& home) const {
    for (auto ni : next)
        if (!ni.assigned())
            return true;
    return false;
}

Choice* FocacciTSPBrancher::choice(const Space&, Archive& e) {
    int index, value;
    e >> index >> value;
    return new FocacciTSPBranchChoice(*this, index, value);
}

ExecStatus FocacciTSPBrancher::commit(Space& home, const Choice& c, unsigned int a) {
    const FocacciTSPBranchChoice& bc = static_cast<const FocacciTSPBranchChoice&>(c);
    const int index = bc.next_index;
    const int value = bc.next_value;

    if (a == 0)
        return me_failed(next[index].eq(home,value)) ? ES_FAILED : ES_OK;
    else
        return me_failed(next[index].nq(home,value)) ? ES_FAILED : ES_OK;

}

void FocacciTSPBrancher::print(
    const Space& home,
    const Choice& c,
    unsigned int a,
    ostream& o) const {

    const FocacciTSPBranchChoice& bc = static_cast<const FocacciTSPBranchChoice&>(c);
    int index = bc.next_index;
    int value = bc.next_value;

    if (a == 0)
        o << "(" << problem.nodes[index] << " " << problem.nodes[value] << ")";
}

unsigned int FocacciTSPBrancher::closest_feasible_arc_index(unsigned int from, unsigned int index) {
    while (index < problem.arcs_size(from)) {
        if (next[from].in((int) problem.arc(from, index).to_index))
            break;
        ++index;
    }
    return index;
}
