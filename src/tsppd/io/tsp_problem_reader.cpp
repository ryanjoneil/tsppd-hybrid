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

#include <regex>
#include <string>

#include <boost/algorithm/string.hpp>

#include <tsppd/io/tsp_problem_reader.h>
#include <tsppd/util/exception.h>
#include <tsppd/util/string.h>

using namespace TSPPD::Data;
using namespace TSPPD::IO;
using namespace TSPPD::Util;
using namespace std;

TSPPDProblem TSPProblemReader::read(string filename) {
    // Defaults for problem attributes.
    string name = "unknown";
    auto type = TSP;
    string comment = "no comment";
    unsigned int dimension = 0;
    auto edge_weight_type = EUC_2D;
    vector<string> nodes{};
    vector<pair<double, double>> coordinates{};

    // Read in the input problem definition.
    ifstream is(filename);
    if (!is.is_open())
        throw TSPPDException("unable to read file '" + filename + "'");

    vector<string> lines;
    string line;
    while (getline(is, line)) {
        boost::trim_left(line);
        boost::trim_right(line);
        lines.push_back(line);
    }

    bool in_edge_weight = false;
    bool in_coordinate = false;
    bool in_precedence = false;

    vector<string> edge_weight_lines;
    vector<string> coordinate_lines;
    vector<string> precedence_lines;

    for (auto line : lines) {
        if (boost::starts_with(line, "EOF")) {
            break;

        } else if (boost::starts_with(line, "EDGE_WEIGHT_SECTION")) {
            in_edge_weight = true;
            in_coordinate = false;
            in_precedence = false;

        } else if (boost::starts_with(line, "NODE_COORD_SECTION")) {
            in_edge_weight = false;
            in_coordinate = true;
            in_precedence = false;

        } else if (boost::starts_with(line, "PRECEDENCE_SECTION")) {
            in_edge_weight = false;
            in_coordinate = false;
            in_precedence = true;

        } else if (in_edge_weight) {
            edge_weight_lines.push_back(line);

        } else if (in_coordinate) {
            coordinate_lines.push_back(line);

        } else if (in_precedence) {
            precedence_lines.push_back(line);

        } else {
            // Parse problem metadata.
            auto comps = resplit(line, regex("\\s*:\\s*"));

            if (comps.size() != 2)
                throw TSPPDException("invalid input '" + line + "'");

            auto lhs = comps[0];
            auto rhs = comps[1];

            if (lhs == "NAME") {
                name = rhs;
            } else if (lhs == "COMMENT") {
                comment = rhs;
            } else if (lhs == "TYPE") {
                if (rhs != "TSP") throw TSPPDException("invalid type '" + rhs + "'");
            } else if (lhs == "DIMENSION") {
                dimension = stoi(rhs);
            } else if (lhs == "EDGE_WEIGHT_TYPE") {
                if (rhs != "EUC_2D" && rhs != "EXPLICIT") throw TSPPDException("invalid edge weight type '" + rhs + "'");
                if (rhs == "EXPLICIT") edge_weight_type = EXPLICIT;
            } else if (lhs == "EDGE_WEIGHT_FORMAT" && edge_weight_type == EXPLICIT) {
                if (rhs != "LOWER_DIAG_ROW") throw TSPPDException("invalid edge weight format '" + rhs + "'");
            } else {
                throw TSPPDException("invalid field '" + lhs + "'");
            }
        }
    }

    is.close();

    // Parse explicit edge weights
    auto edge_weights = read_edge_weights(edge_weight_lines);

    // Parse node and coordinate section.
    auto coord_data = read_coordinate(coordinate_lines);
    nodes = coord_data.first;
    coordinates = coord_data.second;

    // Pickup and delivery pairs, if they exist.
    auto precedence = read_precedence(precedence_lines);

    // Make sure the problem can actually be solved.
    TSPPDProblem problem(name, type, comment, dimension, edge_weight_type, nodes, coordinates, precedence, edge_weights);
    problem.validate();

    return problem;
}

vector<vector<int>> TSPProblemReader::read_edge_weights(vector<string> lines) {
    vector<vector<int>> weights;
    for (auto line : lines) {
        vector<int> w;
        for (auto weight : wssplit(line))
            w.push_back(stoi(weight));
        weights.push_back(w);
    }
    return weights;
}

pair<vector<string>, vector<pair<double, double>>> TSPProblemReader::read_coordinate(vector<string> lines) {
    vector<string> nodes;
    vector<pair<double, double>> coordinates;

    for (auto line : lines) {
        boost::trim_left(line);
        boost::trim_right(line);
        auto comps = wssplit(line);
        if (comps.size() != 3)
            throw TSPPDException("invalid coordinate input '" + line + "'");

        nodes.push_back(comps[0]);
        coordinates.push_back({stod(comps[1]), stod(comps[2])});
    }

    return {nodes, coordinates};
}

vector<pair<string, string>> TSPProblemReader::read_precedence(vector<string> lines) {
    vector<pair<string, string>> pairs;

    for (auto line : lines) {
        boost::trim_left(line);
        boost::trim_right(line);
        auto comps = wssplit(line);
        if (comps.size() != 2)
            throw TSPPDException("invalid precedence input '" + line + "'");

        pairs.push_back({comps[0], comps[1]});
    }

    return pairs;
}
