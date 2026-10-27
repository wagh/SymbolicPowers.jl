# Worked examples

These files are independent, runnable Julia scripts for mathematicians using
SymbolicPowers.jl. Each script constructs its own rings and ideals, explains
the mathematical example, prints results, and checks expected answers with
`@assert`. They use public OSCAR and SymbolicPowers APIs only.

From a clone of the package, instantiate its environment once:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Run one example, or all examples in a single Julia process:

```sh
julia --project=. examples/quickstart.jl
julia --project=. examples/curves.jl
julia --project=. examples/run_all.jl
```

The commands above assume the package root is the working directory. In a
Julia REPL started with `julia --project=.`, one can instead use
`include("examples/containment.jl")`. An absolute path also works with
`include`. The all-examples runner isolates each script in a separate module;
no results from an earlier script are needed by another.

## Choose a topic

Start with `quickstart.jl`. The remaining files can be read in any order.

| File | Mathematical examples | Main functions |
|:--|:--|:--|
| [quickstart.jl](quickstart.jl) | Coordinate points, an embedded prime, a finite field | A first tour of powers, containment and invariants |
| [symbolic_power_methods.jl](symbolic_power_methods.jl) | Seven characteristic-zero methods on the triangle ideal; a complete-intersection point | `symbolic_power`, `is_symbolic_equal_ordinary` |
| [embedded_and_mixed_ideals.jl](embedded_and_mixed_ideals.jl) | Embedded components; incomparable minimal primes of different heights; nonhomogeneous ideals | `associated_prime_heights`, `big_height`, `minimal_part`, both symbolic conventions |
| [projective_points.jl](projective_points.jl) | Four rational points, a double point, an embedded irrelevant component | `symbolic_power` with automatic and saturation methods |
| [curves.jl](curves.jl) | The curve `(t^3,t^4,t^5)` with an explicit element of `P^(2) \ P^2`; the twisted cubic | Primary and unmixed powers, `initial_degree`, `min_degree_symbolic_power` |
| [joins.jl](joins.jl) | Joins of axes; a radical ideal; why nonradical joins differ from symbolic powers | `symbolic_join`, `symbolic_power_join` |
| [positive_characteristic.jl](positive_characteristic.jl) | Frobenius versus ordinary powers, complete intersections, a monomial curve, characteristic two | `frobenius_power`, `symbolic_power_frobenius`, `symbolic_power_join` |
| [containment.jl](containment.jl) | Least symbolic exponent, largest ordinary exponent, a minimal-prime answer of zero | `containment_problem`, `is_symbolic_power_contained_in_power`, `is_symbolic_equal_ordinary` |
| [polyhedra_and_waldschmidt.jl](polyhedra_and_waldschmidt.jl) | Triangle vertices, nonsquarefree lattice holes, exact constants and sampled bounds | `symbolic_polyhedron`, `waldschmidt_constant`, `waldschmidt_upper_bound`, `lower_bound_resurgence` |
| [defect_and_regularity.jl](defect_and_regularity.jl) | A missing symbolic generator, linear complete intersections, nonconstant regularity ratios | `symbolic_defect`, `symbolic_regularity_samples`, `asymptotic_regularity` |
| [packing.jl](packing.jl) | A packed four-cycle, a triangle, a Konig graph with a non-Konig minor | `squarefree_gens`, `squarefree_in_codim`, `is_konig`, `is_packed`, `non_packed_substitutions` |
| [fields_and_gradings.jl](fields_and_gradings.jl) | `QQ(sqrt(2))`, `GF(9)`, positive weights, a multigrading | Powers over different fields, Frobenius coefficients, weighted degree invariants |
| [edge_cases.jl](edge_cases.jl) | Zero/unit ideals, zeroth powers, no variables, rejected hypotheses | Boundary conventions and informative `ArgumentError`s |

Together the examples exercise every exported function. The package CI runs
`run_all.jl`, so the mathematical assertions are also checked after changes.

## Interpret the results

- The default convention uses all associated primes. With
  `use_minimal_primes=true`, even the first symbolic power can differ from
  the input ideal. The embedded-ideal and containment examples show why.
- Equality of ideals does not require the same displayed generators. Use
  `I == J`, rather than equality of their generator arrays. Generator order
  and printed representations can vary between OSCAR versions.
- `waldschmidt_constant` is exact for monomial ideals.
  `waldschmidt_upper_bound` and `lower_bound_resurgence` return certified
  bounds; a finite sample alone does not establish an exact value.
- `symbolic_regularity_samples` returns finite data.
  `asymptotic_regularity(...; approximate=true)` returns their minimum ratio,
  which is not a certified asymptotic invariant or a general bound on it.
- Join and Frobenius helpers check their radicality/primality hypotheses.
  Skip a check only when that hypothesis has been established independently.
- Degree invariants require a positive integer grading. Symbolic powers can
  also be computed for nonhomogeneous ideals and in multigraded rings;
  regularity samples require the standard grading.

The examples use small powers so the computations stay manageable. A first
run also loads OSCAR and compiles Julia methods; it can take substantially
longer than later computations. These scripts are tutorials, not timing
benchmarks; matched Macaulay2 benchmarks live in
[`benchmark/comparison`](../benchmark/comparison).
