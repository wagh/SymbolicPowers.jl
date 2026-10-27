```@meta
CurrentModule = SymbolicPowers
```

# Macaulay2 comparison

Measured 26 October 2026 (UTC) on the unmodified GitHub release **0.3.4**, commit
`7cbd6f6daee81b589868cdd74bed0bbd058307df`. Julia 1.10.12 / OSCAR 1.8.2 was
compared with Macaulay2 1.26.06 / SymbolicPowers 2.0 on one Linux cloud host.
These timings describe that release, not the subsequent review fixes in 0.3.5.

The final review candidate 0.3.5 passes 19,875 tests and agrees with all 58
M2 benchmark answers. Its certified complete-intersection shortcut reduces
`frobenius_linear101_sp2` from 950.71 ms to 3.8902 ms on the same host;
the native M2 median was 282.03 ms. These are separate candidate measurements.
The audit data archive retains all five candidate samples for every group.

## Procedure

Each operation has one untimed warm-up, then five fresh-input calls. The
tables report median wall time in milliseconds. Ring/ideal construction,
Macaulay2 `MonomialIdeal` conversion, package loading and Julia compilation
are excluded. The systems ran sequentially with `JULIA_NUM_THREADS=1`,
`OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`; M2 used `--no-threads` and
Julia used `--startup-file=no --pkgimages=no`.

There are 40 ordinary-input groups and 18 additional specialized-input
groups. Every answer agreed: ideal generators were transferred as exact
coefficients/exponents and compared by ideal equality in OSCAR, outside the
timers; scalar answers were compared exactly. Julia uses `MPolyIdeal`; the
second table uses `MonomialIdeal` in M2. Those conversions occur before timing,
so these are distinct input-representation experiments. All specialized
groups are standard graded. M2's ordinary rings carry the standard grading;
Julia's undecorated rings are interpreted as standard graded for degree
invariants. The weighted curve explicitly has weights `(3,4,5)` in both.

The exact constructors, operations, result checks and reproduction commands
are in `benchmark/comparison` in the checkout. Every original sample and
minimum/maximum is retained in the audit data archive. Medians close to one
another are not reliable universal rankings. In particular, M2 samples show
garbage-collection/system variability in some groups. No aggregate speedup
over unrelated cases is asserted. The release's historical 128-case totals
cannot be reconstructed from the original checkout, which did not contain
that benchmark's inputs and harness; this is a separate reproducible suite.

## Ordinary M2 inputs

| Case | M2 (ms) | Julia 0.3.4 (ms) | Lower median |
|:--|--:|--:|:--|
| `triangle_sp6` | 1.0521 | 0.13413 | Julia |
| `triangle_sp12` | 1.079 | 0.10222 | Julia |
| `fermat_sp2` | 18.6 | 1.298 | Julia |
| `fermat_sp3` | 19.248 | 3.0598 | Julia |
| `curve_sp3` | 197.96 | 6.8377 | Julia |
| `curve_weighted_sp3` | 15.542 | 0.97748 | Julia |
| `twisted_sp3` | 269.39 | 4.892 | Julia |
| `embedded_sp3` | 1.8157 | 0.083956 | Julia |
| `embedded_min_sp3` | 12.669 | 0.13586 | Julia |
| `ci_monomial_sp8` | 1.2504 | 0.060552 | Julia |
| `affine_sp3` | 53.249 | 8.9155 | Julia |
| `join_triangle4` | 17.694 | 19.457 | M2 |
| `join_triangle4_core` | 18.9 | 1.1726 | Julia |
| `join_twisted3_core` | 28.363 | 3.3398 | Julia |
| `frobenius_curve3` | 254.7 | 10.945 | Julia |
| `frobenius_curve4` | 263.24 | 24.575 | Julia |
| `contain_triangle8` | 144.8 | 0.29353 | Julia |
| `contain_fermat2` | 371.46 | 10.62 | Julia |
| `equal_triangle4` | 13.948 | 0.062033 | Julia |
| `height_fermat3` | 273.05 | 6.4942 | Julia |
| `wald_triangle` | 149.21 | 14.102 | Julia |
| `wald_cycle5` | 313.93 | 16.548 | Julia |
| `defect_triangle2` | 1.9722 | 0.058988 | Julia |
| `resurgence_triangle5` | 225.29 | 15.347 | Julia |
| `reg_triangle4` | 4.8888 | 4.4643 | Julia |
| `konig_complete6` | 148.91 | 0.007992 | Julia |
| `konig_hyper6` | 897.16 | 0.027902 | Julia |
| `packed_cycle4` | 269.14 | 0.10561 | Julia |
| `packed_cycle5` | 11.097 | 0.014202 | Julia |
| `packed_triangle_ambient8` | 992.77 | 0.012168 | Julia |
| `squarefree_m5` | 8.9829 | 0.056636 | Julia |
| `squarefree_m6` | 205.4 | 0.081262 | Julia |
| `points8_sp2` | 25.654 | 5.1551 | Julia |
| `points8_sp4` | 463.15 | 33.33 | Julia |
| `fatpoints4_sp2` | 14.649 | 1.1406 | Julia |
| `points8_reg3` | 218.99 | 25.019 | Julia |
| `fermat_ci3` | 998.73 | 30.517 | Julia |
| `mixed_sp2` | 389.7 | 26.717 | Julia |
| `frobenius_linear101_sp2` | 282.03 | 950.71 | M2 |
| `frobenius_point5_sp2` | 52.652 | 6.5851 | Julia |

## Specialized M2 inputs

Julia still receives a fresh `MPolyIdeal`. The `mono_` prefix identifies the
M2 input representation, not a different ideal or operation. Repeated Julia
measurements differ slightly because they are separate timed groups.

| Case | M2 (ms) | Julia 0.3.4 (ms) | Lower median |
|:--|--:|--:|:--|
| `mono_triangle_sp6` | 0.65186 | 0.057086 | Julia |
| `mono_triangle_sp12` | 0.94724 | 0.069204 | Julia |
| `mono_embedded_sp3` | 1.3491 | 0.065689 | Julia |
| `mono_ci_monomial_sp8` | 0.82865 | 0.037887 | Julia |
| `mono_contain_triangle8` | 22.939 | 0.21579 | Julia |
| `mono_equal_triangle4` | 8.9345 | 0.049765 | Julia |
| `mono_wald_triangle` | 159.35 | 11.675 | Julia |
| `mono_wald_cycle5` | 328.28 | 12.847 | Julia |
| `mono_defect_triangle2` | 2.8756 | 0.048793 | Julia |
| `mono_resurgence_triangle5` | 228.12 | 13.609 | Julia |
| `mono_reg_triangle4` | 4.2424 | 3.5672 | Julia |
| `mono_konig_complete6` | 12.503 | 0.007381 | Julia |
| `mono_konig_hyper6` | 68.804 | 0.03337 | Julia |
| `mono_packed_cycle4` | 280.9 | 0.060882 | Julia |
| `mono_packed_cycle5` | 12.753 | 0.008322 | Julia |
| `mono_packed_triangle_ambient8` | 1252.4 | 0.007371 | Julia |
| `mono_squarefree_m5` | 2.9822 | 0.059931 | Julia |
| `mono_squarefree_m6` | 11.805 | 0.05955 | Julia |

## Input and interpretation details

- `triangle`: `(xy,xz,yz)` over `QQ`.
- `fermat`: `(x(y^3-z^3),y(z^3-x^3),z(x^3-y^3))`; `fermat_ci3`
  explicitly requests the component complete-intersection method in both.
- `curve`: `(y^2-xz,x^2*y-z^2,x^3-yz)`, the curve `(t^3,t^4,t^5)`.
  The weighted version uses `(3,4,5)`, and Frobenius curve cases use `GF(3)`.
- `twisted`: the twisted cubic in four variables.
- `embedded`: `(x^3,x^2*y^2,x*z^3,y^4,y^2*z)`; `_min` uses only minimal primes.
- `ci_monomial`: `(x^2,y^2)`; `affine`: `(bc+1,abd)` in four variables.
- `points8`: the points `(1:0:0),(0:1:0),(0:0:1),(1:1:1),(1:2:3),
  (3:1:-2),(2:3:1),(1:-1:2)` over `QQ`. `fatpoints4` doubles the first
  of the first four points and leaves the other three reduced.
- `mixed`: `(a^2+b^2-c) ∩ (a-1,b-c^2,d)` over `QQ`.
- Frobenius linear example: `(x+y,z)` over `GF(101)` at order two.
  Frobenius point example: `(x-1,y,z+x)` over `GF(5)` at order two.
- `contain_*8` asks for the least symbolic exponent contained in the eighth
  ordinary power. `equal_*4` compares symbolic/ordinary fourth powers.
- `resurgence_*5` is the lower bound with five samples and the Waldschmidt
  option. `reg_*` compares the minimum of the indicated number of finite
  regularity ratios, not an exact asymptotic invariant.
- `konig_hyper6` is the squarefree cubic Veronese ideal on six variables;
  `konig_complete6` is the ideal of six variables. `packed_cycle4/5` uses
  graph edge ideals; `packed_triangle_ambient8` has five inactive variables.
  `squarefree_m5/6` uses the ideal of five/six variables.

Julia's public radicality/primality checks are included except in `_core`
rows. M2's join helper has no radicality check. Thus `join_triangle4` includes
different validation work; the `_core` row is the fair comparison of the join
computations with established hypotheses.

## Correctness and scope differences

M2 supports some quotient ambient rings; this Julia package accepts polynomial
ideals over fields. For example, M2 computed `(a,c)^(2)=(a)` in
`QQ[a,b,c]/(ab-c^2)` in the correctness controls. Both packages expose the
all-associated-prime and minimal-prime conventions.

M2 agreement is not the only correctness test. Separate native M2 controls
reproduce these discrepancies, which are excluded from performance rankings:

- For `I=(x^3,xy)` in `QQ[x,y]`, the minimal-prime symbolic square is `(x^2)`,
  which is not contained in `I`. M2's containment predicate nevertheless
  returns `true`, applying the all-associated-prime height bound to the wrong
  convention. Its reverse containment search also returns 1 where the correct
  largest ordinary exponent is 0.
- For `J=((x+y)y,(x+y)z,yz)` in four variables and `I=J ∩ m^3`, the default
  symbolic square is `I^2`: the irrelevant maximal ideal already belongs to
  `Ass(R/I)`. M2's general fallback drops a new height-three component of
  `I^2`, because it tests equality with old associated primes instead of
  containment in them. Julia keeps it.
- M2's `symbolicPower(I,0)` does not return the unit ideal in the control;
  Julia explicitly defines the zeroth power to be the unit ideal.
- For weights `(2,3)` and `I=(x^3,y^2)`, M2's `Ideal` Waldschmidt overload
  returns 6 but the `MonomialIdeal` overload returns 2 for the same ideal.
  Julia consistently returns the weighted value 6. Specialized M2 timings
  above therefore do not include weighted invariants.
- M2's unchecked join helper is documented for prime input; applying it to
  `(x^2)` produces `(x^3)` at order two, not the symbolic square `(x^4)`.
  Julia rejects nonradical input by default. This is a hypothesis boundary,
  not an in-domain counterexample to M2's documented prime helper.

The Julia API also distinguishes exact monomial Waldschmidt constants,
certified sampled upper bounds, and finite regularity data. A finite minimum
of regularity ratios is not generally a bound on asymptotic regularity.

## Cases where M2 was faster

1. **Frobenius power of `(x+y,z)` over `GF(101)`, order two**:
   M2 **282.03 ms**, Julia 0.3.4 **950.71 ms**, about **3.37×** in M2's favor.
   Both validate primality. At height two the colon exponent is 199; M2 uses
   one quotient by `P^199`, whereas 0.3.4 uses 199 successive quotients.
   Version 0.3.5 adds a sufficient complete-intersection certificate and
   returns `P^2` directly for this input, with primality validation retained.
2. **Triangle join at order four, with Julia's default check**:
   the recorded M2 median is **17.69 ms**, versus **19.46 ms** in Julia.
   This modest difference includes extra validation only on Julia's side and
   substantial variation in M2's samples. With `check=false` and radicality
   independently known, Julia's join core takes **1.17 ms**, versus M2's
   **18.90 ms** in the separate core group. It is not an M2 algorithm advantage.

No other group had a lower M2 median in this suite. These results do not
establish that Julia is faster for every other input.
