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

#ifndef TSPPD_UTIL_STRING_H
#define TSPPD_UTIL_STRING_H

#include <regex>
#include <string>
#include <vector>

namespace TSPPD {
    namespace Util {
        // resplit(string, regex) splits a string by a regular expression
        std::vector<std::string> resplit(std::string s, std::regex re) {
            std::vector<std::string> comps{};
            std::regex_token_iterator<std::string::iterator> iter(s.begin(), s.end(), re, -1);
            std::regex_token_iterator<std::string::iterator> end;
            while (iter != end) {
                comps.push_back(*iter);
                ++iter;
            }
            return comps;
        }

        // resplit(string, regex) splits a string by whitespace
        std::vector<std::string> wssplit(std::string s) {
            return resplit(s, std::regex("[\\s\\t]+"));
        }
    }
}

#endif
