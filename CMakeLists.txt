cmake_minimum_required(VERSION 3.7.2)
project(tsppd)

option(GIST "GIST" ON)
if(${GIST})
    add_definitions(-DGIST)
endif(${GIST})

# Speed up builds with pre-compiled headers
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")
include(cotire)

# Compiler flags
set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -std=c++14 -Werror -Wall -Wextra -Wno-c++98-compat -Wno-unused-parameter -pedantic-errors")

# Where to look for header filess
include_directories(src/)
if(APPLE)
    include_directories(/usr/local/include/)
    include_directories(/Library/gurobi800/mac64/include/)
endif(APPLE)

# Where to look for libraries
if(APPLE)
    link_directories(/usr/local/lib/)
    link_directories(/Library/gurobi800/mac64/lib/)
endif(APPLE)

# Libraries we link against
set(BOOST_LIBRARIES boost_program_options)
if(${GIST})
    set(GECODE_LIBRARIES gecodefloat gecodegist gecodeint gecodekernel gecodeminimodel gecodesearch gecodeset gecodesupport)
else(${GIST})
    set(GECODE_LIBRARIES gecodefloat gecodeint gecodekernel gecodeminimodel gecodesearch gecodeset gecodesupport)
endif(${GIST})
set(GUROBI_LIBRARIES gurobi_c++ gurobi80 m)
set(ALL_LIBRARIES ${BOOST_LIBRARIES} ${GECODE_LIBRARIES} ${GUROBI_LIBRARIES})

# This allows "make install" to put all the headers in the right place.
set(TSPPD_LIB_HEADERS
    src/tsppd/ap/ap_solver.h
    src/tsppd/ap/gurobi_ap_solver.h
    src/tsppd/ap/primal_dual_ap_solver.h
    src/tsppd/data/tsppd_arc.h
    src/tsppd/data/tsppd_problem.h
    src/tsppd/data/tsppd_problem_generator.h
    src/tsppd/data/tsppd_search_statistics.h
    src/tsppd/data/tsppd_solution.h
    src/tsppd/data/tsppd_tree.h
    src/tsppd/io/tsp_problem_reader.h
    src/tsppd/io/tsp_problem_writer.h
    src/tsppd/io/tsp_solution_writer.h
    src/tsppd/solver/enumerative/enumerative_tsp_solver.h
    src/tsppd/solver/enumerative/enumerative_tsppd_solver.h
    src/tsppd/solver/focacci/brancher/focacci_tsp_branch_choice.h
    src/tsppd/solver/focacci/brancher/focacci_tsp_brancher.cpp
    src/tsppd/solver/focacci/brancher/focacci_tsp_closest_neighbor_brancher.h
    src/tsppd/solver/focacci/brancher/focacci_tsp_regret_brancher.h
    src/tsppd/solver/focacci/brancher/focacci_tsp_sequential_closest_neighbor_brancher.h
    src/tsppd/solver/focacci/dual/focacci_closest_neighbor_dual.h
    src/tsppd/solver/focacci/dual/focacci_tsp_dual.h
    src/tsppd/solver/focacci/propagator/focacci_tsppd_assignment_propagator.h
    src/tsppd/solver/focacci/propagator/focacci_tsppd_omc_propagator.h
    src/tsppd/solver/focacci/propagator/focacci_tsppd_precede_cost_propagator.h
    src/tsppd/solver/focacci/propagator/focacci_tsppd_precede_propagator.h
    src/tsppd/solver/focacci/propagator/focacci_tsppd_precede_set_propagator.h
    src/tsppd/solver/focacci/focacci_tsp_solver.h
    src/tsppd/solver/focacci/focacci_tsp_space.h
    src/tsppd/solver/focacci/focacci_tsppd_solver.h
    src/tsppd/solver/focacci/focacci_tsppd_space.h
    src/tsppd/solver/oneil/oneil_tsppd_callback.h
    src/tsppd/solver/oneil/oneil_tsppd_solver.h
    src/tsppd/solver/ruland/callback/ruland_precedence_callback.h
    src/tsppd/solver/ruland/callback/ruland_subtour_elimination_callback.h
    src/tsppd/solver/ruland/callback/ruland_tsp_callback.h
    src/tsppd/solver/ruland/callback/ruland_tsp_callback_handler.h
    src/tsppd/solver/ruland/ruland_subtour_finder.h
    src/tsppd/solver/ruland/ruland_tsp_solver.h
    src/tsppd/solver/ruland/ruland_tsppd_plus_solver.h
    src/tsppd/solver/ruland/ruland_tsppd_solver.h
    src/tsppd/solver/sarin/sarin_tsp_callback.h
    src/tsppd/solver/sarin/sarin_tsp_solver.h
    src/tsppd/solver/sarin/sarin_tsppd_solver.h
    src/tsppd/solver/tsp_solver.h
    src/tsppd/util/exception.h
    src/tsppd/util/stacktrace.h
    src/tsppd/util/string.h)

set(TSPPD_LIB_SOURCES
    src/tsppd/ap/ap_solver.cpp
    src/tsppd/ap/gurobi_ap_solver.cpp
    src/tsppd/ap/primal_dual_ap_solver.cpp
    src/tsppd/data/tsppd_problem.cpp
    src/tsppd/data/tsppd_problem_generator.cpp
    src/tsppd/data/tsppd_solution.cpp
    src/tsppd/data/tsppd_tree.cpp
    src/tsppd/io/tsp_problem_reader.cpp
    src/tsppd/io/tsp_problem_writer.cpp
    src/tsppd/io/tsp_solution_writer.cpp
    src/tsppd/solver/enumerative/enumerative_tsp_solver.cpp
    src/tsppd/solver/enumerative/enumerative_tsppd_solver.cpp
    src/tsppd/solver/focacci/brancher/focacci_tsp_brancher.cpp
    src/tsppd/solver/focacci/brancher/focacci_tsp_closest_neighbor_brancher.cpp
    src/tsppd/solver/focacci/brancher/focacci_tsp_regret_brancher.cpp
    src/tsppd/solver/focacci/brancher/focacci_tsp_sequential_closest_neighbor_brancher.cpp
    src/tsppd/solver/focacci/dual/focacci_closest_neighbor_dual.cpp
    src/tsppd/solver/focacci/propagator/focacci_tsppd_assignment_propagator.cpp
    src/tsppd/solver/focacci/propagator/focacci_tsppd_omc_propagator.cpp
    src/tsppd/solver/focacci/propagator/focacci_tsppd_precede_cost_propagator.cpp
    src/tsppd/solver/focacci/propagator/focacci_tsppd_precede_set_propagator.cpp
    src/tsppd/solver/focacci/focacci_tsp_solver.cpp
    src/tsppd/solver/focacci/focacci_tsp_space.cpp
    src/tsppd/solver/focacci/focacci_tsppd_solver.cpp
    src/tsppd/solver/focacci/focacci_tsppd_space.cpp
    src/tsppd/solver/oneil/oneil_tsppd_callback.cpp
    src/tsppd/solver/oneil/oneil_tsppd_solver.cpp
    src/tsppd/solver/ruland/callback/ruland_precedence_callback.cpp
    src/tsppd/solver/ruland/callback/ruland_subtour_elimination_callback.cpp
    src/tsppd/solver/ruland/callback/ruland_tsp_callback_handler.cpp
    src/tsppd/solver/ruland/ruland_subtour_finder.cpp
    src/tsppd/solver/ruland/ruland_tsp_solver.cpp
    src/tsppd/solver/ruland/ruland_tsppd_plus_solver.cpp
    src/tsppd/solver/ruland/ruland_tsppd_solver.cpp
    src/tsppd/solver/sarin/sarin_tsp_callback.cpp
    src/tsppd/solver/sarin/sarin_tsp_solver.cpp
    src/tsppd/solver/sarin/sarin_tsppd_solver.cpp
    src/tsppd/util/stacktrace.cpp)

# tsppd library
add_library(tsppd SHARED ${TSPPD_LIB_SOURCES})
target_link_libraries(tsppd ${ALL_LIBRARIES})
cotire(tsppd)

# ap test executable
set(TSPPD_BIN_SOURCES src/tsppd-bin/ap.cpp)
add_executable(ap-bin ${TSPPD_BIN_SOURCES})
target_link_libraries(ap-bin tsppd ${ALL_LIBRARIES})
set_target_properties(ap-bin PROPERTIES OUTPUT_NAME ap)
cotire(ap-bin)

# tsppd executable
set(TSPPD_BIN_SOURCES src/tsppd-bin/main.cpp)
add_executable(tsppd-bin ${TSPPD_BIN_SOURCES})
target_link_libraries(tsppd-bin tsppd ${ALL_LIBRARIES})
set_target_properties(tsppd-bin PROPERTIES OUTPUT_NAME tsppd)
cotire(tsppd-bin)
