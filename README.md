# SymbolicPowers.jl

Symbolic powers of ideals in polynomial rings over fields, for
[OSCAR](https://www.oscar-system.org).

The package computes symbolic powers by several algorithms, studies the
containment problem (for which m is I^(m) ⊆ I^n?), and computes invariants such as
the Waldschmidt constant, the symbolic defect and bounds for the resurgence.
It also implements the packing problem for squarefree monomial ideals. Exact
results, proven bounds and finite samples are kept apart in the API.

## Installation

The package requires Julia 1.10 or later and OSCAR 1.8 or later (OSCAR runs on
Linux and macOS, and on Windows through WSL). From a local clone:

```julia
using Pkg
Pkg.develop(path = "/path/to/SymbolicPowers.jl")
```

or, from the published repository:

```julia
Pkg.add(url = "https://github.com/wagh/SymbolicPowers.jl")
```

## Quick start

```julia
using Oscar, SymbolicPowers

R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
I = ideal(R, [x*y, x*z, y*z])                 # three coordinate points of P^2

symbolic_power(I, 2) == I^2 + ideal(R, [x*y*z])   # true
containment_problem(I, 2)                     # 3: least m with I^(m) ⊆ I^2
waldschmidt_constant(I)                       # 3//2 (exact)
waldschmidt_upper_bound(I; sample_size = 1)   # 2   (an upper bound)
lower_bound_resurgence(I; use_waldschmidt = true)  # 4//3
symbolic_defect(I, 2)                         # 1
```

## Definition

By default (Hochster–Huneke)

    I^(n) = ⋂_{P ∈ Ass(R/I)} (I^n R_P ∩ R).

With `use_minimal_primes = true` the intersection runs over the minimal primes
of `I` only. The two definitions agree when `I` has no embedded primes (for
instance when `I` is radical), but can differ even for `n = 1`:

```julia
R, (x, y) = polynomial_ring(QQ, [:x, :y])
I = ideal(R, [x^2, x*y])                                        # (x) ∩ (x^2, y)
symbolic_power(I, 2) == I^2                                     # true
symbolic_power(I, 2; use_minimal_primes = true) == ideal(R, [x^2])  # true
```

The zeroth symbolic power is the unit ideal.

## Algorithms

`symbolic_power(I, n; method = ...)`:

| `method` | Requirements (checked) | Algorithm |
|:--|:--|:--|
| `:auto` | – | monomial engine, dimension-one shortcuts, unmixed hull, or componentwise localization |
| `:primary` | – | primary decomposition of `I^n`, keep the components localizing correctly |
| `:monomial` | monomial ideal | powers of the localized components (Cooper–Embree–Hà–Hoefel) |
| `:saturation` | homogeneous, `dim R/I = 1` | complete-intersection or certified rational fat-point shortcut; otherwise saturate at the irrelevant ideal, respecting embedded primes |
| `:unmixed` | all associated primes of equal height | equidimensional hull of `I^n` |
| `:complete_intersection` | components are complete intersections of height `codim(I)` | intersect their powers |
| `:join` | radical; characteristic 0 or finite field | Zariski–Nagata: the join `I ⋆ m^n` |
| `:frobenius` | prime; characteristic `p > 0` | `P^[q] : P^(h(q-1)-n+1)` with `q = p^e ≥ n` |

A method whose requirements fail throws an `ArgumentError`; no other method is
substituted silently.

Both symbolic conventions use the monomial and unmixed shortcuts. Under the
minimal-prime convention, automatic dimension-one selection uses `:saturation`
only when the complete-intersection or rational fat-point shortcut applies;
otherwise it uses the componentwise method. This retains isolated components
of every height, unlike taking one equidimensional hull of a mixed ideal.

## Exact values, bounds and samples

- `waldschmidt_constant` is exact and requires a monomial ideal (it minimizes
  the degree over the symbolic polyhedron). For other homogeneous ideals,
  `waldschmidt_upper_bound` returns `min α(I^(n))/n` over a finite range, which
  is a proven upper bound.
- `lower_bound_resurgence` returns a proven lower bound for the resurgence.
- `symbolic_regularity_samples` returns `reg(I^(n))` for finitely many `n`;
  such samples do not bound the asymptotic regularity in general.
  `asymptotic_regularity` (the minimum ratio) requires `approximate = true`.
- `symbolic_polyhedron`: exponents of monomials in `I^(n)` lie in `n` times the
  polyhedron; the converse holds for squarefree monomial ideals only.

## Functions

| Area | Functions |
|:--|:--|
| Symbolic powers | `symbolic_power`, `symbolic_join`, `symbolic_power_join`, `symbolic_power_frobenius`, `frobenius_power` |
| Invariants | `big_height`, `associated_prime_heights`, `minimal_part`, `initial_degree`, `min_degree_symbolic_power` |
| Containment | `is_symbolic_equal_ordinary`, `is_symbolic_power_contained_in_power`, `containment_problem`, `lower_bound_resurgence` |
| Asymptotics | `symbolic_polyhedron`, `waldschmidt_constant`, `waldschmidt_upper_bound`, `symbolic_defect`, `symbolic_regularity_samples`, `asymptotic_regularity` |
| Packing | `squarefree_gens`, `squarefree_in_codim`, `is_konig`, `is_packed`, `non_packed_substitutions` |

Every function has a docstring (`?symbolic_power` in the REPL). The manual in
`docs/src` explains the mathematics and proves the correctness of each
algorithm.

## Limitations

- Ideals must live in polynomial rings over fields; quotient rings, local
  rings and power series rings are not supported.
- Degree invariants need homogeneous ideals in a positive `Z`-grading;
  regularity needs the standard grading.
- Primary decomposition and elimination can be expensive, and OSCAR's primary
  decomposition can fail or fail to terminate over very small prime fields. There, use `:join`
  or `:frobenius` (with `check = false` if radicality or primality is known).
- Over number fields, redundant backend primary components are removed before
  their radicals are used as associated primes; normalization is cached.
- Nontrivial ideal powers inherit the exponent and resource limits of the
  OSCAR/Singular backend. Fat-point exponent products are checked for integer
  overflow. The monomial engine uses widened integers for degree comparisons;
  degree invariants returning `Int` still raise an overflow error if their
  result does not fit.
- `is_packed` examines up to `3^a` minors, where `a` is the number of variables
  occurring in the generators.

## Tests and documentation

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The tests compare all methods with each other, check monomial symbolic powers
against an independent exponent-vector oracle, compare the combinatorial
packing code and the containment searches with brute-force computations, and
cover the zero and unit ideals, embedded primes, weighted and multigraded
rings, number fields (including redundant primary decompositions), independent
Hasse-derivative tests, and positive characteristic. The CI also checks
JuliaFormatter 1.0.62 formatting with OSCAR's settings.

See [the Macaulay2 comparison](docs/src/comparison.md) for measured examples,
the complete benchmark procedure and the cases where Macaulay2 was faster.
Benchmark scripts are in [`benchmark/comparison`](benchmark/comparison).
Times depend on the coefficient field, generators, input representation,
chosen method and whether the ideal object is reused.

## References

- M. Hochster, C. Huneke, *Comparison of symbolic and ordinary powers of ideals*, Invent. Math. 147 (2002).
- S. Sullivant, *Combinatorial symbolic powers*, J. Algebra 319 (2008).
- S. Cooper, R. Embree, H. T. Hà, A. Hoefel, *Symbolic powers of monomial ideals*, Proc. Edinb. Math. Soc. 60 (2017).
- H. Dao, A. De Stefani, E. Grifo, C. Huneke, L. Núñez-Betancourt, *Symbolic powers of ideals*, arXiv:1708.03010.
- B. Drabkin, E. Grifo, A. Seceleanu, B. Stone, *Calculations involving symbolic powers*, J. Softw. Algebra Geom. 9 (2019).

The package is inspired by the Singular library `symbolicPowers.lib` and the
Macaulay2 package `SymbolicPowers`. License: GPL-3.0-or-later.
