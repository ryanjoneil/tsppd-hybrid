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

#include <iostream>

#include <execinfo.h>
#include <signal.h>
#include <unistd.h>

#include <tsppd/util/stacktrace.h>

namespace TSPPD {
    namespace Util {
        void dump_stack(int sig) {
            fprintf(stderr, "Error: signal %d:\n", sig);

            void *array[20];
            size_t size;

            // get void*'s for all entries on the stack
            size = backtrace(array, sizeof(array));

            // print out all the frames to stderr
            backtrace_symbols_fd(array, size, STDERR_FILENO);
        }

        void stacktrace_dumper(int sig) {
            dump_stack(sig);
            exit(1);
        }

        void enable_stacktraces() {
            signal(SIGSEGV, stacktrace_dumper);
            signal(SIGBUS, stacktrace_dumper);
            signal(SIGINT, stacktrace_dumper);
            signal(SIGQUIT, stacktrace_dumper);
            signal(SIGILL, stacktrace_dumper);
            signal(SIGABRT, stacktrace_dumper);
            signal(SIGFPE, stacktrace_dumper);
            signal(SIGTERM, stacktrace_dumper);
            signal(SIGSYS, stacktrace_dumper);

            signal(SIGUSR1, dump_stack);
        }
    }
}
