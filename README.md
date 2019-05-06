# TSPPD Hybrid Optimization Code

This source code accompanies the following papers:

* Exact Methods for Solving Traveling Salesman Problems with Pickup and Delivery in Real Time
* Integer Models for the Traveling Salesman Problem with Pickup and Delivery
* Decision Diagrams for Solving Traveling Salesman Problems with Pickup and Delivery in Real Time

The `src/` directory contains C++14 source code for the various TSPPD models
used in the paper. All source is licensed under the ZIB Academic License, and
can only be used or referenced for academic purposes. There is no warranty and
your mileage may vary. See the file `tsppd/LICENSE`.

Building `tsppd` requires `cmake`, `clang` or `g++`, and the following
libraries:

* Boost >= 1.6.7
* Gecode 6.2.0
* Gurobi 8.1.1

In order to build it, create a directory called build and use `cmake`.

```
mkdir build
cd build
cmake ..
make
```

If you don't have GIST support in your copy of Gecode, append "-DGIST=off" to
the cmake invocation.

In the build succeeds, you should have a tsppd binary in your tsppd/build/
directory. You can run it against the problem instances using the `-i` flag
to specify the input data and `-s` to specify the solver.

```
./tsppd -i <input file> -s <solver>
```

For example:

```
./tsppd -i file.tsp -s tsppd-enum
```

Help is available by using the `-h` flag.

```
tsppd command line interface:
  -h [ --help ]              produce help message
  -H [ --no-header ]         do not print csv header line
  -a [ --asymmetric ]        asymmetric mode (.tsp only) -
                             upper triangular matrix = lower * U(0.7,1.3)
  -s [ --solver ] arg        solver slug
  -i [ --input ] arg         input tsplib file
  -f [ --format ] arg        output format: {human|csv} (default=human)
  -r [ --random-seed ] arg   random seed (default=0)
  -n [ --random-size ] arg   randomly generated instance with n pairs
  -w [ --save ] arg          save problem instance to (modifed) tsplib format
                             file
  -o [ --solver-option ] arg solver option (e.g. foo=bar)
  -p [ --threads ] arg       threads (default=1)
  -t [ --time-limit ] arg    time limit in millis
  -l [ --soln-limit ] arg    stop after n solutions are found
```

Solver-specific options are passed using the -o flag. Multiple options can be
specified by using more than one -o flag (e.g. -o foo=bar -o baz=qux).
Solver-specific options follow. Not all of these are used in the papers.
Solvers denoted `atsppd-*` use asymmetric integer formulations.

```
atsppd-ap
    sec:      subtour elimination constraint type
              - cutset:  x(delta(S)) >= 1
              - subtour: sum { i,j in S } in x_{i,j} <= |S| - 1 (default)

atsppd-oneil
    relax:    relax model and add sec as violated {on|off} (default=off)
    sec:      subtour elimination constraint type for relaxed models
              - cutset:  x(delta(S)) >= 1
              - subtour: sum { i,j in S } in x_{i,j} <= |S| - 1 (default)

atsppd-sarin
    prec:     relaxed prec form that uses x or y variables {x|y} (default=x)
    relax:    relax model, add sec and prec as violated {on|off} (default=off)
    sec:      subtour elimination constraint type for relaxed models
              - cutset:  x(delta(S)) >= 1
              - subtour: sum { i,j in S } in x_{i,j} <= |S| - 1 (default)
              - y:       y_ij + x_ji + y_jk + y_ki <= 2

tsppd-focacci
    brancher: branching scheme {cn, regret, seq-cn} (default=regret)
    dl:       discrepancy limit (lds only)
    filter:   variable domain filtering mechanism (default=none)
              - ap:   assignment problem reduced cost propagator
              - hk:   1-tree bound and marginal cost propagator
              - aphk: additive bounding using ap + hk
              - hkap: additive bounding using hk + ap
    gist:     enables interactive search tool (implies search=bab)
    hk-iter:  max iterations for hk 1-tree bound (default=10)
    omc:      order matching constraints (default=off)
    precede:  precedence propagator type {set, cost, all} (default=set)
    search:   search engine {bab, dfs, lds} (default=bab)

tsppd-ruland
    sec:      subtour elimination constraint type
              - cutset:  x(delta(S)) >= 2
              - subtour: sum { i,j in S } in x_{i,j} <= |S| - 1 (default)
              - hybrid:  subtour if |S| <= (N + 1) / 3, else cutset

tsppd-ruland+, atsppd-oneil+, atsppd-sarin+
    warm-time: milliseconds spent in tsppd-focacci warm start
    warm-soln: solution limit for warm start
```

When warm starting MIP using the CP solver, as in `tsppd-ruland+`, `tsppd-focacci`
solver options are passed on to the CP solver.
