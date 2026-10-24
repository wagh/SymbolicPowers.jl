```@meta
CurrentModule = SymbolicPowers
DocTestSetup = quote
  using Oscar, SymbolicPowers
end
```

# SymbolicPowers.jl

`SymbolicPowers` computes symbolic powers of ideals in polynomial rings over
fields with [OSCAR](https://www.oscar-system.org), together with the
containment problem, Waldschmidt constants, resurgence bounds, symbolic
defects, symbolic polyhedra of monomial ideals and the packing problem.

## Installation

The package needs Julia 1.10 or later and OSCAR 1.8 or later. From a local
clone:

```julia
using Pkg
Pkg.develop(path = "/path/to/SymbolicPowers.jl")
using Oscar, SymbolicPowers
```

## A first computation

```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal(R, [x*y, x*z, y*z]);   # three coordinate points of P^2

julia> symbolic_power(I, 2) == I^2 + ideal(R, [x*y*z])
true

julia> containment_problem(I, 2)
3

julia> waldschmidt_constant(I)
3//2

julia> lower_bound_resurgence(I; use_waldschmidt = true)
4//3

julia> symbolic_defect(I, 2)
1
```

## Scope

- Ideals in polynomial rings over fields (for example `QQ`, number fields and
  finite fields), graded or not. Quotient rings, localizations, power series
  rings and coefficients in `ZZ` are not supported.
- Degree invariants (initial degree, Waldschmidt constants, symbolic defect)
  require homogeneous ideals in a positive `Z`-grading; regularity samples
  require the standard grading.
- Exact results, proven bounds and finite samples have different functions:
  `waldschmidt_constant` (exact, monomial ideals) versus
  `waldschmidt_upper_bound` (a bound), and `symbolic_regularity_samples`
  (finite data, no asymptotic claim).

The [Mathematical guide](@ref) states the definitions and conventions, and
[Algorithms and guarantees](@ref) proves the correctness of every method.
The changes between versions are listed in `CHANGELOG.md`.

## Package overview

```@docs
SymbolicPowers
```
