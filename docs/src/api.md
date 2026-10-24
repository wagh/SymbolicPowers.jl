```@meta
CurrentModule = SymbolicPowers
DocTestSetup = quote
  using Oscar, SymbolicPowers
end
```

# API

## Powers and joins

```@docs
symbolic_power
symbolic_join
symbolic_power_join
symbolic_power_frobenius
frobenius_power
```

## Associated primes and degrees

```@docs
big_height
associated_prime_heights
minimal_part
initial_degree
min_degree_symbolic_power
```

## Containment

```@docs
is_symbolic_equal_ordinary
is_symbolic_power_contained_in_power
containment_problem
lower_bound_resurgence
```

## Polyhedra, defect and asymptotic information

```@docs
symbolic_polyhedron
waldschmidt_constant
waldschmidt_upper_bound
symbolic_defect
symbolic_regularity_samples
asymptotic_regularity
```

## Packing

```@docs
squarefree_gens
squarefree_in_codim
is_konig
is_packed
non_packed_substitutions
```
