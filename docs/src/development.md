```@meta
CurrentModule = SymbolicPowers
DocTestSetup = quote
  using Oscar, SymbolicPowers
end
```

# Development and validation

## Running the tests

From the package root:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

The test suite contains

- examples of the Macaulay2 package `SymbolicPowers` and of
  `symbolicPowers.lib`, with applicable methods cross-checked on representative
  families;
- an independent oracle for monomial ideals: associated primes are found as
  monomial colon ideals, localizations are computed on exponent vectors, and
  membership of every monomial in a box is checked for both conventions;
- brute-force references, built only from OSCAR ideal operations, for the
  combinatorial König and packing code and for the containment searches;
- boundary cases: zero and unit ideals, zeroth and first powers, embedded and
  mixed-height primes, nonradical joins, Frobenius exponents, overflow,
  weighted, zero-weight and multigraded rings, nonrational points and finite
  extension fields.
- number fields, including a deliberately redundant but valid decomposition
  supplied through OSCAR's cache to exercise its documented backend contract;
- independent Hasse-derivative membership tests, including characteristic two
  with symbolic exponent larger than the characteristic, and an imperfect
  coefficient field where the join method must be rejected;
- a nonmonomial ideal whose square has a new associated prime contained in an
  old embedded associated prime; selecting primes by equality loses it;
- squarefree products compared with direct ideal powering on small examples,
  and a 20-variable case for which the direct power would be impractical.

## Building the manual

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

This runs all doctests, and `checkdocs = :exports` makes sure that every
exported function appears in the manual.

## Conventions

The code follows the OSCAR style guide: two-space indentation, snake_case
names, `is_`/`has_` prefixes for predicates, `@req` for argument checks (which
throw `ArgumentError`), and `@doc raw"""..."""` docstrings with `jldoctest`
examples. `.JuliaFormatter.toml` contains OSCAR's formatter settings.

The package has been formatted with JuliaFormatter 1.0.62 (the version-1
formatter recommended by OSCAR). To reproduce the formatting check from the
package root:

```sh
julia --startup-file=no -e 'using Pkg; Pkg.activate(; temp=true); Pkg.add(PackageSpec(name="JuliaFormatter", version="1.0.62")); using JuliaFormatter; @assert format(".")'
```

The CI workflow runs this check as well as package tests and the manual build.
If formatting a modified file returns `false`, rerun until no changes remain.

Internal helpers start with an underscore. The keyword `check = false` of
`symbolic_power_join` and `symbolic_power_frobenius` skips a hypothesis test;
it must only be used when the hypothesis is known to hold.

## Performance notes

- Reuse the same ideal object when computing several exponents. The package
  caches minimal and associated primes, the defining primes for each symbolic
  convention, and localized components on that object. The two conventions
  have separate cache keys. These caches supplement OSCAR's decomposition
  caches without modifying them. Number-field decompositions are normalized
  once and cached separately.
- Structural caches are independent of the exponent. The package does not
  retain every computed symbolic power indefinitely. As with OSCAR's own
  cached ideal invariants, do not mutate an ideal's generators or its cached
  mathematical data after computing invariants; construct a new ideal instead.
- Containment searches gallop from `n` and cache symbolic powers, so the
  largest symbolic power computed is at most about twice the answer;
  `lower_bound_resurgence` shares the cache between different `n`.
- `is_symbolic_equal_ordinary` checks direct containment for monomial ideals,
  avoiding primary decomposition of the ordinary power. Other ideals use the
  associated-prime criterion.
- Monomial symbolic defects count minimal symbolic generators outside the
  ordinary power, avoiding quotient-module construction. Other ideals use a
  minimal generating set in the quotient ring.
- Squarefree symbolic polyhedra are built directly from their coordinate-prime
  inequalities. Waldschmidt constants request an exact rational linear-program
  minimum rather than collecting and inspecting every vertex in Julia.
- `is_konig`, `is_packed` and `non_packed_substitutions` work on the supports
  of the minimal generators; no ideal arithmetic is needed for the minors.
- The initial degree of a homogeneous ideal is read off homogeneous
  generators; a Gröbner basis is only computed if the given generators are not
  homogeneous. A separate unit-ideal test is deliberately avoided here,
  because OSCAR's `is_one` can itself compute a Gröbner basis.
- `squarefree_in_codim` enumerates disjoint products and minimizes their
  supports, instead of constructing all generators of the ordinary power.
- Radical ideals are joined with ``\mathfrak m^n`` directly, without splitting
  into minimal primes.

`benchmark/smoke.jl` prints a few timings.

Measure both fresh ideals and repeated calls on the same ideal after an
untimed warm-up. Julia compilation and package loading can dominate a first
call, while fresh ideals do not benefit from structural caches. Report these
costs separately from steady-state algebra. Proper positive weights can also
make a weighted-homogeneous problem substantially easier for the backend.
