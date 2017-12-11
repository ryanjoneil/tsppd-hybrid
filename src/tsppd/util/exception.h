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

#ifndef TSPPD_UTIL_EXCEPTION_H
#define TSPPD_UTIL_EXCEPTION_H

#include <exception>
#include <string>

namespace TSPPD {
    namespace Util {
        class TSPPDException : public std::exception {
        public:
            TSPPDException(const std::string message) : std::exception(), message(message) { }
            const char* what() const noexcept {
                return message.c_str();
            }

        private:
            const std::string message;
        };
    }
}

#endif
