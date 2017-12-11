#!/bin/bash

set -e

PATH=$PATH:"$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/build

if [ "$PROB" = "" ]; then
    PROB="tsp tsppd"
fi

if [ "$SEED" = "" ]; then
    ITER_SEQ=$(seq 0 9)
    SEED=$RANDOM
    UPDATE_SEED=1
else
    ITER_SEQ=0
    UPDATE_SEED=0
fi

if [ "$SIZE" = "" ]; then
    ITER_SIZE=$(seq 2 6)
else
    ITER_SIZE=$SIZE
fi

run() {
    $1 -f csv | csvcut -c primal | tail -1 >> .out
}

FAILED=0
for PROB in $PROB; do
    for SIZE in $ITER_SIZE; do
        for ITER in $ITER_SEQ; do
            rm -f .out

            SEED_OUT=$(printf "%05d" $SEED)
            CMD="tsppd -n $SIZE -r $SEED"

            /bin/echo -n "[$SIZE / $SEED_OUT] testing $PROB solvers... "

            # Enum
            if [ $SIZE -le 5 ]; then
                run "$CMD -s $PROB-enum"
            fi

            # CP
            run "$CMD -s $PROB-cp"
            run "$CMD -s $PROB-cp -o threads=2"

            for BRANCH in cn seq-cn regret; do
                run "$CMD -s $PROB-cp -o brancher=$BRANCH"
            done

            if [ "$PROB" == "tsppd" ]; then
                run "$CMD -s $PROB-cp -o ap=on"

                for PRECEDE in all cost set; do
                    run "$CMD -s $PROB-cp -o precede=$PRECEDE"
                done
            fi

            # MIP
            run "$CMD -s $PROB-mip"
            for SEC in cutset subtour hybrid; do
                run "$CMD -s $PROB-mip -o sec=$SEC"
            done
            if [ "$PROB" == "tsppd" ]; then
                run "$CMD -s $PROB-mip+ -o warm-time=50"
            fi

            if [ $(sort .out | uniq | wc -l) -eq 1 ]; then
                echo "ok"
            else
                echo "failed"
                FAILED=1
            fi

            if [ $UPDATE_SEED -eq 1 ]; then
                SEED=$RANDOM
            fi
        done
    done
done

if [ $FAILED -eq 1 ]; then
    echo "tests failed"
    exit 1
else
    echo "tests passed"
fi
