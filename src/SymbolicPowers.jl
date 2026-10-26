"""
    SymbolicPowers

Symbolic powers of ideals in polynomial rings over fields, built on OSCAR.

The package is inspired by the Singular library `symbolicPowers.lib`
(V. Wagh) and the Macaulay2 package `SymbolicPowers` (E. Grifo et al.):
several algorithms for computing symbolic powers (primary decomposition,
saturation, unmixed part, joins, Frobenius powers, monomial and
complete-intersection shortcuts), the containment problem, the packing problem
for monomial ideals, symbolic polyhedra, the Waldschmidt constant, resurgence
bounds, symbolic defect, and finite symbolic-regularity samples. Exact monomial
Waldschmidt constants are distinguished from sampled upper bounds.

Unless stated otherwise, the `n`-th symbolic power of `I` is

```math
I^{(n)} = \\bigcap_{P \\in \\operatorname{Ass}(R/I)} \\left(I^n R_P \\cap R\\right)
```

(Hochster–Huneke). Passing `use_minimal_primes = true` switches to the
definition that only uses the minimal primes of `I`. Both definitions agree for
ideals without embedded primes, in particular for radical ideals.
"""
module SymbolicPowers

using Oscar

export associated_prime_heights,
  asymptotic_regularity,
  big_height,
  containment_problem,
  frobenius_power,
  initial_degree,
  is_konig,
  is_packed,
  is_symbolic_equal_ordinary,
  is_symbolic_power_contained_in_power,
  lower_bound_resurgence,
  min_degree_symbolic_power,
  minimal_part,
  non_packed_substitutions,
  squarefree_gens,
  squarefree_in_codim,
  symbolic_defect,
  symbolic_join,
  symbolic_polyhedron,
  symbolic_power,
  symbolic_power_frobenius,
  symbolic_power_join,
  symbolic_regularity_samples,
  waldschmidt_constant,
  waldschmidt_upper_bound

include("utils.jl")
include("monomial.jl")
include("symbolic_power.jl")
include("points.jl")
include("radical.jl")
include("join.jl")
include("containment.jl")
include("packing.jl")
include("asymptotic.jl")

end # module
