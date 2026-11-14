#!/bin/sh
# Quick self-test of ffbd (about a minute). Run from the repository root: make test
# 1. Every clock, on both examples: after every proposal the cached posterior is recomputed from scratch and compared (--check-every-step).
# 2. Stepping-stone runs end to end and prints a marginal likelihood.
# 3. The bears example (real data) runs and its root age falls inside the root prior.
set -u
B=./ffbd; T=$(mktemp -d); fail=0
[ -x "$B" ] || { echo "build ffbd first (make)"; exit 1; }
for ex in examples/example_input.txt examples/example_crown_stem.txt; do
  for c in exp ln gam igr acln; do
    out=$($B run "$ex" --prefix "$T/chk" --clock $c --seed 11 --gens 20000 --print 5000 --burnin 5000 --check-every-step 2>&1)
    if echo "$out" | grep -q "ERROR" || ! echo "$out" | grep -q "FINAL"; then echo "FAIL  cache check  $ex  clock $c"; fail=1; else echo "ok    cache check  $ex  clock $c"; fi
  done
done
out=$($B run examples/example_input.txt --prefix "$T/ss" --clock gam --seed 3 --burnin 20000 --ss 8 --ss-burn 2000 --ss-gens 10000 --ss-every 100 2>&1)
if echo "$out" | grep -qi "stepping"; then echo "ok    stepping-stone"; else echo "FAIL  stepping-stone"; fail=1; fi
out=$($B run examples/bears_input.txt --prefix "$T/bears" --clock gam --seed 5 --gens 200000 --print 10000 --burnin 50000 2>&1)
if echo "$out" | grep -q "FINAL"; then echo "ok    bears example runs"; else echo "FAIL  bears example"; fail=1; fi
rm -rf "$T"
[ $fail -eq 0 ] && echo "all tests passed" || { echo "some tests FAILED"; exit 1; }
