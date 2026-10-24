```@meta
CurrentModule = SymbolicPowers
DocTestSetup = quote
  using Oscar, SymbolicPowers
end
```

# Mathematical guide

## Two definitions

Let ``R=K[x_1,\ldots,x_d]`` and let ``I`` be a proper ideal. The default is

```math
I^{(n)}_{\mathrm{Ass}}=
\bigcap_{P\in\operatorname{Ass}(R/I)}(I^nR_P\cap R).
```

With `use_minimal_primes=true`, use

```math
I^{(n)}_{\mathrm{Min}}=
\bigcap_{P\in\operatorname{Min}(R/I)}(I^nR_P\cap R).
```

They agree if there are no embedded primes, including for radical ideals.
The package sets every zeroth symbolic power to ``R``. For ``n\geq1`` the
zero ideal has symbolic power zero and the unit ideal has symbolic power ``R``.
The all-associated first symbolic power is ``I``; the minimal-prime first
symbolic power can be strictly larger and is `minimal_part(I)`.

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> I = ideal(R, [x^2, x*y]);

julia> symbolic_power(I, 2) == I^2
true

julia> symbolic_power(I, 1; use_minimal_primes=true) == ideal(R, [x])
true

julia> symbolic_power(I, 2; use_minimal_primes=true) == ideal(R, [x^2])
true
```

`minimal_part` keeps **all** isolated primary components, even if the minimal
primes have different heights. It differs from an equidimensional hull, which
retains only components of maximal dimension.

## Containment

`is_symbolic_power_contained_in_power(I,m,n)` tests
``I^{(m)}\subseteq I^n`` and respects either definition. The all-associated
case uses the uniform containment theorem ``I^{(hn)}\subseteq I^n``, where
``h=\operatorname{bight}(I)``. That shortcut is not applied directly to the
minimal-prime definition of a nonradical ideal.

`containment_problem(I,n)` returns the least ``m`` with that containment.
With `in_symbolic=true`, it returns the greatest ``k`` such that
``I^{(n)}\subseteq I^k``. This can be zero for the minimal-prime definition.
The ordinary direction obeys `max_power` including at the initial exponent.
With `in_symbolic=true`, the search lies in ``0\leq k\leq n`` and
`max_power` does not change that search.

## Degrees and Waldschmidt constants

Degree routines require homogeneous ideals in a positive integer grading.
An undecorated polynomial ring uses the standard grading. For weights
``w_i>0``, the degree of ``x^a`` is ``\sum_i w_i a_i``.
The zero ideal has no finite initial degree and raises `ArgumentError`.

For a fixed symbolic-power convention the filtration is multiplicative:
``I^{(a)}I^{(b)}\subseteq I^{(a+b)}``. Hence initial degrees are subadditive,
so

```math
\widehat\alpha(I)=\lim_{n\to\infty}\frac{\alpha(I^{(n)})}{n}
=\inf_{n\geq1}\frac{\alpha(I^{(n)})}{n}.
```

Every finite minimum of these ratios is an **upper bound**, not a proof of the
exact value. `waldschmidt_upper_bound` computes this bound. Exact monomial
values are obtained by `waldschmidt_constant` using the symbolic polyhedron.

```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal(R, [x*y, x*z, y*z]);

julia> (waldschmidt_upper_bound(I; sample_size=1), waldschmidt_constant(I))
(2, 3//2)
```

## Symbolic polyhedra

For a monomial ideal and a maximal associated prime ``P``, let
``Q_{\subseteq P}=IR_P\cap R``. Its symbolic polyhedron is
``\mathcal S(I)=\bigcap_P\operatorname{NP}(Q_{\subseteq P})``, with
``\operatorname{NP}(Q)=\operatorname{conv}(\text{generator exponents})+
\mathbb R_{\geq0}^d``.
Then the monomial exponent set satisfies
``\mathcal L(I^{(n)})\subseteq n\mathcal S(I)\cap\mathbb N^d``.
Equality is valid for squarefree monomial ideals but is false in general.
For example, ``(1,1)\in\operatorname{NP}((x^2,y^2))`` while
``xy\notin(x^2,y^2)``. This distinction concerns integral closure and cannot
be repaired by merely increasing numerical precision.

The minimum of a positive degree functional on ``\mathcal S(I)`` is still
exactly the monomial Waldschmidt constant. A rational optimum has a positive
multiple lying in the exponent sets of every local component power, giving
attainment along an infinite subsequence.

With the minimal-prime flag the polyhedron uses the corresponding isolated
components and the same argument applies.

## Symbolic defect and regularity

`symbolic_defect(I,n)` is the minimal number of homogeneous module generators
of ``I^{(n)}/I^n``. It is computed as an ideal in ``R/I^n``. It is not the
difference between the numbers of generators of the two ideals.

`symbolic_regularity_samples` uses standard Castelnuovo-Mumford regularity
of the ideal as a module, not the regularity of its quotient. It requires the
standard grading. A finite sample minimum is not generally a bound on
``\lim_n\operatorname{reg}(I^{(n)})/n``. Such regularity sequences need not
be subadditive. The API makes no general assertion that a limit exists.

## Packing

Nonconstant monomials form a regular sequence exactly when their supports are pairwise
disjoint. `is_konig` searches for such a family of size ``\operatorname{ht}(I)``.
For squarefree ideals this is the matching-versus-cover property of the
corresponding clutter. `is_packed` checks that every minor has that property.
The Conforti-Cornuejols assertion relating packing to equality of all symbolic
and ordinary powers is a conjecture and is not used as an algorithm.
