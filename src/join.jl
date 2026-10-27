###############################################################################
#  Joins of ideals (port of `symbolicJoin` from symbolicPowers.lib)
###############################################################################

@doc raw"""
    symbolic_join(I::MPolyIdeal, J::MPolyIdeal) -> MPolyIdeal

Return the join ``I \star J`` of two ideals of ``R = K[x_1, \dots, x_d]``,
computed by eliminating the variables ``t`` from
``I(t) + J(x - t) \subseteq K[x, t]``.

This is defined for arbitrary ideals. Its relation to symbolic powers
(``I \star \mathfrak m^n = I^{(n)}``) requires `I` to be radical; see
[`symbolic_power_join`](@ref).

See Definition 2.16 in Dao, De Stefani, Grifo, Huneke and Núñez-Betancourt,
*Symbolic powers of ideals* (arXiv:1708.03010), and Sullivant, *Combinatorial
symbolic powers*, J. Algebra 319 (2008).

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> symbolic_join(ideal(R, [x, y]), ideal(R, [x, z])) == ideal(R, [x])
true
```
"""
function symbolic_join(I::MPolyIdeal, J::MPolyIdeal)
  R = _check_ideal(I)
  @req base_ring(J) === R "the ideals must live in the same ring"
  d = ngens(R)
  d == 0 && return I + J # Over K itself there are no variables to eliminate.
  K = coefficient_ring(R)
  names = vcat([Symbol("_sp_x", i) for i in 1:d], [Symbol("_sp_t", i) for i in 1:d])
  S, v = polynomial_ring(K, names; cached=false)
  x = v[1:d]
  t = v[(d + 1):(2d)]
  in_t = hom(R, S, t)                    # f(x) -> f(t)
  in_x_minus_t = hom(R, S, x .- t)       # f(x) -> f(x - t)
  big = ideal(S, vcat([in_t(f) for f in gens(I)], [in_x_minus_t(f) for f in gens(J)]))
  E = eliminate(big, t)
  back = hom(S, R, vcat(gens(R), [zero(R) for _ in 1:d]))
  return ideal(R, [back(f) for f in gens(E)])
end

@doc raw"""
    symbolic_power_join(I::MPolyIdeal, n::Int; check::Bool = true) -> MPolyIdeal

Return ``I \star \mathfrak m^n`` where ``\mathfrak m = (x_1, \dots, x_d)``.
For a radical ideal `I` of a polynomial ring over a perfect field this is the
symbolic power ``I^{(n)}``: by the Zariski–Nagata theorem ``I^{(n)}`` consists
of the polynomials all of whose (Hasse) derivatives of order less than `n` lie
in `I`, and this is exactly the elimination ideal computed by the join. This is
the method of `symbolicPowerofPrimeIdeal` in `symbolicPowers.lib`, applied to
the radical ideal directly instead of to each of its minimal primes.

The coefficient field must have characteristic zero or be an OSCAR finite
field (`FinField`), sufficient conditions for perfectness. With
`check = true` (default) radicality of `I` is
verified; use `check = false` only if `I` is known to be radical. For a
non-radical ideal the join is not the symbolic power: for `I = (x^2)` in
``\mathbb Q[x, y]`` and `n = 2` the join is `(x^3)`, whereas ``I^{(2)} = (x^4)``.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal(R, [x*y, x*z, y*z]);

julia> symbolic_power_join(I, 2) == I^2 + ideal(R, [x*y*z])
true
```
"""
function symbolic_power_join(I::MPolyIdeal, n::Int; check::Bool=true)
  R = _check_ideal(I)
  _check_power(n)
  K = coefficient_ring(R)
  @req characteristic(K) == 0 || K isa Oscar.FinField "the join method requires a perfect field (characteristic 0 or an OSCAR finite field)"
  # A polynomial ring in zero variables is a field: its only ideals are 0
  # and R, both radical. The elimination backend requires an indeterminate.
  ngens(R) == 0 && return n == 0 ? _unit_ideal(R) : I
  check && @req is_radical(I) "the join method requires a radical ideal"
  n == 0 && return _unit_ideal(R)
  n == 1 && return I
  (is_zero(I) || _is_unit_ideal(I)) && return I
  return symbolic_join(I, _maximal_homogeneous_ideal(R)^n)
end

###############################################################################
#  Frobenius method for prime ideals in positive characteristic
###############################################################################

@doc raw"""
    symbolic_power_frobenius(P::MPolyIdeal, n::Int; check::Bool = true) -> MPolyIdeal

Return the symbolic power ``P^{(n)}`` of a prime ideal `P` in a polynomial ring
over a field of characteristic ``p > 0``, as the colon ideal
``P^{[q]} : P^{h(q-1) - n + 1}``, where ``h`` is the height of `P` and ``q`` is
the smallest power of ``p`` with ``q \geq n``.

Since ``P^{[q]}`` is ``P``-primary (Frobenius is flat on a regular ring), the
colon is ``P``-primary, and in the regular local ring ``R_P`` it equals
``P^n R_P`` whenever ``n \leq q``.

Set `check = false` to skip the primality test.

If the supplied nonzero generators number at most the height of `P`, they
certify that the prime is a complete intersection. Then every ordinary power
is `P`-primary, so `P^n` is returned directly instead of computing the colon.

# Examples
```jldoctest
julia> F, (x, y, z) = polynomial_ring(GF(3), [:x, :y, :z]);

julia> P = ideal(F, [y^2 - x*z, x^2*y - z^2, x^3 - y*z]);   # (t^3, t^4, t^5)

julia> symbolic_power_frobenius(P, 2) == symbolic_power_join(P, 2)
true
```
"""
function symbolic_power_frobenius(P::MPolyIdeal, n::Int; check::Bool=true)
  R = _check_ideal(P)
  _check_power(n)
  p = characteristic(coefficient_ring(R))
  @req p > 0 "method :frobenius requires positive characteristic"
  if ngens(R) == 0
    @req !check || is_zero(P) "method :frobenius requires a prime ideal"
    return n == 0 ? _unit_ideal(R) : P
  end
  check && @req is_prime(P) "method :frobenius requires a prime ideal"
  n == 0 && return _unit_ideal(R)
  (n == 1 || is_zero(P)) && return P
  h = codim(P)
  q = ZZ(p)
  while q < n
    q *= p
  end
  e = ZZ(h) * (q - 1) - n + 1
  @req q <= typemax(Int) && e <= typemax(Int) "the Frobenius exponent is too large for this implementation"
  # A prime generated by height-many elements in a regular polynomial ring
  # is a complete intersection, whose powers remain primary to that prime.
  # Keep the established exponent validation and the primality check above.
  count(!iszero, gens(P)) <= h && return P^n
  Q = frobenius_power(P, Int(q))
  # Compute P^[q] : P^e as ((P^[q] : P) : P) ... : P, using (I : J) : K = I : JK.
  # The power P^e has very many (and redundant) generators, and one quotient
  # by it is much slower than e quotients by P: for (t^3, t^4, t^5) over GF(5)
  # and n = 6 (q = 25, e = 43) 250 s against 5.4 s. Quotients by blocks P^k,
  # k = 2, 3, 4, 6, were slower than k = 1 in every tested case.
  for _ in 1:Int(e)
    Q = quotient(Q, P)
  end
  return Q
end
