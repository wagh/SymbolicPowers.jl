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

- OSCAR caches primary decompositions on ideal objects; helpers such as
  `_associated_primes` and `_defining_primes` rely on this. Number-field
  decompositions are normalized once and cached separately.
- Containment searches gallop from `n` and cache symbolic powers, so the
  largest symbolic power computed is at most about twice the answer;
  `lower_bound_resurgence` shares the cache between different `n`.
- `is_symbolic_equal_ordinary` compares associated primes only.
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
