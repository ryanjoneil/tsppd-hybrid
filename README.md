# TSPPD Solver Code

This source code accompanies the paper "Exact Methods for Solving Traveling
Salesman Problems with Pickup and Delivery in Real Time".

The `src/` directory contains C++14 source code for the various TSPPD models
used in the paper. All source is licensed under the ZIB Academic License, and
can only be used or referenced for academic purposes. There is no warranty and
your mileage may vary. See the file `tsppd/LICENSE`.

Building `tsppd` requires `cmake`, `clang` or `g++`, and the following
libraries:

* Boost >= 1.6.0
* Gecode 6.0.1
* Gurobi 8.0.0

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
directory. You can run it against the problem instances using the -i flag
to specify the input data and -s to specify the solver.

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
  -s [ --solver ] arg        {tsp|tsppd}-{cp|enum|mip}
  -i [ --input ] arg         input tsplib file
  -f [ --format ] arg        output format: {human|csv} (default=human)
  -r [ --random-seed ] arg   random seed (default=0)
  -n [ --random-size ] arg   randomly generated instance with n pairs
  -w [ --save ] arg          save problem instance to (modifed) tsplib
                             format  file
  -o [ --solver-option ] arg solver option (e.g. foo=bar)
  -t [ --time-limit ] arg    time limit in millis
  -l [ --soln-limit ] arg    stop after n solutions are found
```

Solver-specific options are passed using the -o flag. Multiple options can be
specified by using more than one -o flag (e.g. -o foo=bar -o baz=qux).
Solver-specific options used in the paper follow.

```
tsppd-cp solver options:
    brancher: branching scheme {cn, regret, seq-cn} (default=regret)
    precede:  precedence propagator type {set, cost, all} (default=set)
    ap:       assignment problem reduced cost propagator {on, off}
              (default=off)
    threads:  number of threads to use in Gecode (default=1)

tsppd-mip solver options:
    sec: subtour elimination constraint type
        - cutset:      x(delta(S)) >= 2
        - subtour:     sum { i,j in S } in x_{i,j} <= |S| - 1 (default)
        - hybrid:      subtour if |S| <= (N + 1) / 3, else cutset

tsppd-mip+ solver options:
    warm-time: milliseconds spent in tsppd-cp warm start
```

When warm starting MIP using the CP solver, as in `tsppd-mip+`, `tsppd-cp`
solver options are passed on to the CP solver.
