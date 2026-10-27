# Changelog

## 0.3.5

- Expand `examples` with twelve independent topic scripts covering every
  exported function, a mathematical guide and an isolated all-examples runner.
  Check the examples in CI on each supported test configuration.
- Check the product of a symbolic exponent and a fat-point multiplicity before
  powering. Integer wraparound could previously return a smaller, incorrect
  power instead of reporting an exponent overflow.
- Use widened integers for monomial degree sums, divisibility-mask threshold
  products and coordinate-prime degree tests. Individually valid `Int`
  exponents may have a total degree exceeding `typemax(Int)`; wrapping could
  retain redundant generators and misclassify a squarefree ideal.
- Decide the unit-ideal test directly for monomial generators over fields,
  avoiding unnecessary Gröbner bases and their exponent restrictions.
- Add large-exponent regressions and declare the `Random` test dependency.
- Use an overflow-safe upper midpoint in reverse containment searches. Compute
  powers of principal monomial ideals by checked scalar multiplication, so
  sparse boundary cases do not require one multiplication per exponent.
- Handle joins and checked power helpers directly in polynomial rings with
  zero variables, where the only ideals are zero and the unit ideal and the
  elimination backend cannot be used.
- In the Frobenius helper, return the ordinary power when the supplied
  generators certify that the prime is a complete intersection. Primality is
  still checked. This avoids a large number of unnecessary quotients, for
  example 199 of them for a linear prime over `GF(101)` at symbolic order two.
- Update installation and algorithm documentation, explain arithmetic limits,
  and add a reproducible comparison with Macaulay2, including input types and
  cases where Macaulay2 is faster.

## 0.3.4

Performance release. Results are unchanged (generator sets may differ). On a
128-case benchmark shared with the Macaulay2 package `SymbolicPowers` 2.0, the
total time drops from 506 s (0.3.3) to 14 s (Macaulay2: 290 s); no case is
slower than 0.3.3 or than Macaulay2, and all results agree with Macaulay2.

- Monomial ideals: a combinatorial engine on exponent vectors
  (`src/monomial.jl`) replaces Gröbner bases and primary decompositions in the
  monomial paths: symbolic powers, minimal and associated primes, big height,
  heights of associated primes, `minimal_part`, containment searches,
  `is_symbolic_equal_ordinary` and the symbolic defect. 13-24x on the larger
  benchmark cases. The squarefree prime caches no longer share one vector.
- Points (homogeneous, dimension one): over QQ and finite fields, schemes of
  (fat) rational points are recognized by projections, certified, and their
  symbolic powers computed as intersections of powers of the point ideals
  instead of saturating `I^n` (8 points in P^2, n = 6: 25 s -> 0.4 s).
- Method selection: unmixedness is decided by comparing `I` with its
  equidimensional hull, without a primary decomposition of `I`; `big_height`
  and `associated_prime_heights` use this as well. Ideals that are not
  unmixed (and dimension-one ideals under the minimal-prime convention without
  a points shortcut) use a componentwise method: the power of each isolated
  primary component is hulled separately, and only the part below embedded
  defining primes needs a primary decomposition of its power (non-homogeneous
  mixed ideal, n = 5: 55 s -> 0.06 s). For the minimal-prime convention,
  `:auto` now also uses `:saturation` (through its shortcuts) and `:unmixed`.
- `symbolic_power_frobenius` computes `P^[q] : P^e` as `e` successive
  quotients by `P` (curve `(t^3,t^4,t^5)` over GF(5), n = 6: 245 s -> 5 s).
- Intersections of several non-monomial ideals are folded pairwise
  (complete-intersection method, Fermat quintic, n = 4: 22 s -> 1.3 s).
- `symbolic_regularity_samples` passes a minimal generating set to
  `cm_regularity` (8 points in P^2, 4 samples: 82 s -> 0.14 s).
- New regression tests compare every shortcut with the general methods
  (19,756 tests in total).

## 0.3.3

- Cache exponent-independent minimal primes, associated primes, defining
  primes, and localized components on each input ideal, with separate keys
  for the two symbolic conventions and without modifying OSCAR's caches.
- Decide monomial symbolic/ordinary equality by direct containment, avoiding
  primary decomposition of the ordinary power.
- Compute monomial symbolic defects by surviving minimal monomial generators
  instead of constructing a quotient module.
- Assemble squarefree symbolic polyhedra directly from coordinate-prime
  inequalities, and compute Waldschmidt constants with an exact rational LP.
- Add reuse, convention-isolation, quotient-module-reference, polyhedron,
  finite-field, and weighted-grading regression checks. Explain the shortcuts
  and cache assumptions in the manual.

## 0.3.2

- The primary method no longer normalizes the decomposition of `I^n` over
  number fields. Its component selection is valid for every primary
  decomposition, redundant or not; only the defining primes of `I` need an
  irredundant decomposition, and those are still normalized. This avoids a
  quadratic number of intersections and containment tests on the (typically
  large) decomposition of the power.
- New test: redundant components injected into the decomposition of the power
  do not change the result, for both conventions.

## 0.3.1

### Correctness and robustness

- Normalize OSCAR's potentially redundant number-field primary decompositions
  before reading their associated primes or applying the component shortcuts.
  Merge components with equal radicals, remove redundant components, cache the
  result separately, and preserve the backend's original cache.
- Return zero and unit ideals directly for positive symbolic powers, avoiding
  unnecessary backend exponent limits; do the same in the direct join helper.
- Reject unsupported positive-characteristic fields in the join method with
  `ArgumentError`, without calling an unimplemented generic `is_finite` method.
- Handle zero and unit inputs in `squarefree_in_codim` explicitly: the zero
  ideal's codimension is zero, while OSCAR assigns infinite codimension to the
  unit ideal. Both return `[one(R)]`, by the documented boundary conventions.
- Correct the containment-search regression described as having an embedded
  maximal prime: intersection with the square of the maximal ideal left the
  quadratic ideal unchanged; its cube creates the intended example.

### Performance and standards

- Avoid the implicit Gröbner basis caused by testing `is_one` before computing
  initial degree from homogeneous generators.
- Enumerate disjoint squarefree products in `squarefree_in_codim` rather than
  constructing the full ideal power, and remove redundant resulting products.
- Run JuliaFormatter 1.0.62 with OSCAR's configuration on all Julia files and
  add a formatting CI job.
- Use a strict distinctness check with the normal ordering in the brute-force
  matching reference, instead of passing a non-strict comparator to `issorted`.

### Tests and documentation

- Add number-field computations and a controlled redundant-decomposition
  fixture for OSCAR's documented backend behavior.
- Add independent Hasse-derivative checks for joins, an imperfect-field
  counterexample, a new-associated-prime regression, and squarefree-product
  comparisons and a 20-variable example.
- Explain number-field normalization, backend exponent limits, the disjoint
  product algorithm, and how to reproduce the formatting check in the manual.

## 0.3.0

### Correctness and robustness

- `symbolic_power(...; method = :saturation)` and the automatic method choice
  decide homogeneity for the ideal (via a reduced Gröbner basis if necessary),
  not for the given generators.
- The complete-intersection method checks that the primary decomposition has
  one component per associated prime.
- The docstring of `symbolic_power` describes the join method as implemented
  (the radical ideal itself is joined with a power of the maximal ideal).

### API

- `is_symbolic_equal_ordinary` accepts `use_minimal_primes`.
- `non_packed_substitutions` only substitutes variables occurring in a minimal
  generator; witnesses involving other variables were redundant copies.

### Performance

- `containment_problem` uses a galloping search from `n` with cached symbolic
  powers; plain bisection of `[n, hn]` could compute symbolic powers of order
  `hn/2` even when the answer is close to `n`. For `in_symbolic = true` the
  search starts at the Hochster–Huneke bound `floor(n/h)`.
- `lower_bound_resurgence` shares symbolic powers between its searches and
  with the sampled Waldschmidt bound.
- `is_symbolic_equal_ordinary` compares associated primes and computes no
  intersection.
- `is_konig`, `is_packed` and `non_packed_substitutions` are computed
  combinatorially on supports (transversal number = codimension), without
  constructing the minors as ideals.

### Code and documentation

- Argument checks use OSCAR's `@req`.
- Shared helpers for the defining primes and the localized components of
  monomial ideals.
- New tests: brute-force references for the König/packing code and for the
  containment searches, both conventions of `is_symbolic_equal_ordinary`,
  homogeneous ideals with nonhomogeneous generators, a third value of the
  corona regularity example.
- README and manual rewritten for users; review artifacts and generated HTML
  removed from the package.

## 0.2.0

### Correctness

- Correct the zero-ideal containment predicate at symbolic exponent zero.
- Validate radicality in the direct symbolic join helper.
- Restrict Frobenius exponents to characteristic powers, avoiding
  generator-dependent output, and protect exponent arithmetic from overflow.
- Compute Konig's property using disjoint supports, including nonsquarefree
  monomials. Require squarefree input for clutter packing predicates.
- Enforce positive integer grading for degree invariants and standard grading
  for Castelnuovo-Mumford regularity samples.
- Enforce the initial containment search cap and validate sample sizes.
- Correct the lattice-point assertion for general symbolic polyhedra.

### API changes

- `waldschmidt_constant` now returns exact monomial values only; sampling is
  available as `waldschmidt_upper_bound`.
- Add `symbolic_regularity_samples`. The legacy `asymptotic_regularity`
  requires `approximate=true` and explicitly reports only a finite statistic.
- Add the minimal-prime flag to monomial polyhedra and Waldschmidt constants.

### Performance and documentation

- Minimize supplied monomial generators without a Groebner basis.
- Avoid a Groebner computation for initial degree of a homogeneous presentation.
- Compute radical joins directly, without splitting into minimal primes.
- Use binary containment searches with a reusable fixed ideal.
- Enumerate substitutions lazily and skip inactive variables in packing tests.
- Add a Documenter manual, mathematical algorithm justifications, regression
  tests, an independent exponent oracle, examples, and reproducible benchmarks.

## 0.1.0

Original supplied package.
