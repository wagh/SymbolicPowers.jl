# Reproduce the Macaulay2 comparison

Run these commands from the package root. The recorded review uses
Macaulay2 1.26.06 / SymbolicPowers 2.0 and Julia 1.10.12 / OSCAR 1.8.2.
The published measurements describe the unmodified 0.3.4 commit
`7cbd6f6daee81b589868cdd74bed0bbd058307df`; running against another checkout
measures that checkout instead. Pin the desired Julia dependency versions
before measuring. The audit data archive supplies the exact manifests.

```sh
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
mkdir -p benchmark/comparison/results
M2 --silent --stop -q --no-readline --no-threads benchmark/comparison/bench.m2 > benchmark/comparison/results/m2.log 2>&1
julia --startup-file=no --pkgimages=no --project=. benchmark/comparison/bench.jl benchmark/comparison/results/m2.log > benchmark/comparison/results/julia.log 2>&1
M2 --silent --stop -q --no-readline --no-threads benchmark/comparison/correctness.m2
julia --startup-file=no --pkgimages=no --project=. benchmark/comparison/correctness.jl
M2 --silent --stop -q --no-readline --no-threads benchmark/comparison/weighted.m2
```

For the exact original release, create a separate worktree at the commit above,
instantiate its pinned environment, and pass its directory to `--project`
while running these harness files from the current package root. Both scripts
accept `BENCH_CASE=triangle_sp12` (or another case identifier) to run one
group. `BENCH_DEFINITIONS_ONLY=1` loads definitions without running the suite.

There are 40 ordinary-input groups and 18 additional `mono_` groups using
Macaulay2's `MonomialIdeal` representation. Julia always receives a fresh
`MPolyIdeal`. All specialized groups use the standard grading: the Macaulay2
Waldschmidt overloads disagree about weighted degree, so substituting input
types for weighted invariants would not give an equivalent comparison.

Input and ring construction, including `MonomialIdeal` conversion, is outside
the timer. One call warms each operation before five fresh-input calls; the
median wall time is reported, together with every raw sample. The scripts run
in sequence. They include public hypothesis checks, except rows explicitly
labelled `_core`, where Julia's `check=false` is used. Macaulay2's join helper
does not check radicality, so `join_triangle4` compares different validation
work. All coefficient fields, generators and weights appear in the scripts.

Macaulay2 emits ideal generators as coefficient/exponent records. Julia
reconstructs them and compares ideals, independently of generator order or
normalization, outside the timers. Boolean and rational answers are compared
exactly. Missing or different answers make the Julia script fail. Regularity
rows compare the same minimum of finitely many ratios, not a certified
asymptotic value. Known Macaulay2 boundary/convention errors are tested
separately and excluded from timing rankings.

Five samples and these small representative inputs do not establish general
scaling laws. Close medians can change with hardware, backend version, garbage
collection and input representation. Julia `--pkgimages=no` was used for the
recorded runs; loading and compilation are excluded from operation timings.
