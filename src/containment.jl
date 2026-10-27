###############################################################################
#  Containment problem
###############################################################################

@doc raw"""
    is_symbolic_equal_ordinary(I::MPolyIdeal, n::Int;
                               use_minimal_primes::Bool = false) -> Bool

Return whether ``I^{(n)} = I^n``.

For monomial ideals, the test checks ``I^{(n)} \subseteq I^n`` directly; the
reverse containment always holds. For other ideals, the test compares
associated primes: ``I^{(n)} = I^n`` holds exactly when
every associated prime of ``I^n`` is contained in a maximal associated prime of
`I` (respectively, with `use_minimal_primes = true`, when ``I^n`` has no
embedded primes).

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> is_symbolic_equal_ordinary(ideal(R, [x*y, x*z, y*z]), 2)
false

julia> is_symbolic_equal_ordinary(ideal(R, [x*z, y*z]), 2)
true

julia> I = ideal(R, [x^2, x*y]);   # (x) ∩ (x^2, y)

julia> is_symbolic_equal_ordinary(I, 1), is_symbolic_equal_ordinary(I, 1; use_minimal_primes = true)
(true, false)
```
"""
function is_symbolic_equal_ordinary(I::MPolyIdeal, n::Int; use_minimal_primes::Bool=false)
  _check_ideal(I)
  _check_power(n)
  (n == 0 || _is_unit_ideal(I) || is_zero(I)) && return true
  !use_minimal_primes && n == 1 && return true
  if is_monomial(I)
    # Both minimal generating sets are canonical (src/monomial.jl).
    return _mono_symbolic_power_exponents(I, n, use_minimal_primes) ==
           _mono_power_exponents(I, n)
  end
  # The primes of the components kept in I^(n) form a subset S of Ass(I^n)
  # that is closed under going down; then Ass(I^(n)) = S, so I^(n) = I^n iff
  # S = Ass(I^n).
  targets = _defining_primes(I, use_minimal_primes)
  return all(P -> any(M -> is_subset(P, M), targets), _associated_primes(I^n))
end

@doc raw"""
    is_symbolic_power_contained_in_power(I::MPolyIdeal, m::Int, n::Int;
                                         use_minimal_primes::Bool = false) -> Bool

Return whether ``I^{(m)} \subseteq I^n``.

For the default definition of symbolic powers the theorem of Ein–Lazarsfeld–Smith
and Hochster–Huneke, ``I^{(hn)} \subseteq I^n`` with ``h`` the big height of
`I`, is used to answer without computation when ``m \geq hn``. This shortcut is
not used for the minimal-prime definition, for which it can fail.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal(R, [x*y, x*z, y*z]);

julia> is_symbolic_power_contained_in_power(I, 2, 2), is_symbolic_power_contained_in_power(I, 3, 2)
(false, true)
```
"""
function is_symbolic_power_contained_in_power(I::MPolyIdeal, m::Int, n::Int;
  use_minimal_primes::Bool=false)
  _check_ideal(I)
  _check_power(m)
  _check_power(n)
  n == 0 && return true
  _is_unit_ideal(I) && return true
  is_zero(I) && return m > 0
  # I^m ⊆ I^(m); for a nonzero proper ideal I^m ⊆ I^n forces m ≥ n
  m < n && return false
  if !use_minimal_primes && ZZ(m) >= ZZ(big_height(I)) * n
    return true
  end
  if is_monomial(I)
    return _mono_is_subset(
      _mono_symbolic_power_exponents(I, m, use_minimal_primes), _mono_power_exponents(I, n)
    )
  end
  return is_subset(symbolic_power(I, m; use_minimal_primes), I^n)
end

# Symbolic powers of a fixed ideal, computed at most once per exponent.
function _cached_symbolic_power!(cache::Dict{Int,T}, I::T, m::Int,
  use_minimal_primes::Bool) where {T<:MPolyIdeal}
  return get!(() -> symbolic_power(I, m; use_minimal_primes), cache, m)
end

# The smallest m with I^(m) ⊆ I^n, for a nonzero proper ideal I and n ≥ 1.
# The predicate m ↦ [I^(m) ⊆ I^n] is monotone, false for m < n, and (default
# definition) true for m ≥ h*n. A galloping search from n is used: it tests
# n, n+1, n+3, n+7, ... and then bisects the last gap, so the largest symbolic
# power computed is less than about twice the answer.
function _containment_threshold(I::T, n::Int, use_minimal_primes::Bool, max_power::Int,
  cache::Dict{Int,T}) where {T<:MPolyIdeal}
  @req n <= max_power "max_power must be at least n"
  # Monomial ideals: containment by divisibility of minimal generators.
  contained = if is_monomial(I)
    Qe = _mono_power_exponents(I, n)
    m -> _mono_is_subset(
      _mono_exponents(_cached_symbolic_power!(cache, I, m, use_minimal_primes)), Qe
    )
  else
    Q = I^n
    m -> is_subset(_cached_symbolic_power!(cache, I, m, use_minimal_primes), Q)
  end
  certified = false               # is containment at `cap` known by theory?
  cap = max_power
  if !use_minimal_primes
    bound = ZZ(big_height(I)) * n
    if bound <= max_power
      cap = Int(bound)
      certified = true
    end
  end
  lo = n - 1                      # largest exponent known not to work
  offset = 0
  hi = cap
  while true
    m = offset >= cap - n ? cap : n + offset
    if m == cap && certified
      hi = cap
      break
    end
    if contained(m)
      hi = m
      break
    end
    @req m < cap "no containment found up to max_power = $max_power"
    lo = m
    offset = offset >= div(cap - n, 2) ? cap - n : 2 * offset + 1
  end
  while hi - lo > 1
    mid = lo + div(hi - lo, 2)
    if contained(mid)
      hi = mid
    else
      lo = mid
    end
  end
  return hi
end

@doc raw"""
    containment_problem(I::MPolyIdeal, n::Int; in_symbolic::Bool = false,
                        use_minimal_primes::Bool = false,
                        max_power::Int = typemax(Int)) -> Int

With `in_symbolic = false` (default), return the smallest `m` such that
``I^{(m)} \subseteq I^n``. The answer satisfies ``n \le m \le hn``, where ``h``
is the big height of `I` (Hochster–Huneke). For the minimal-prime definition
there is no such a priori bound, and the search stops with an `ArgumentError`
after `max_power`.

With `in_symbolic = true`, return the largest `k` such that
``I^{(n)} \subseteq I^k``. Then ``\lfloor n/h \rfloor \le k \le n`` for the
default definition and ``0 \le k \le n`` for the minimal-prime definition;
`max_power` is not used.

The ideal must be nonzero and proper, and `n` must be positive.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal(R, [x*y, x*z, y*z]);

julia> containment_problem(I, 2)
3

julia> containment_problem(I, 3; in_symbolic = true)
2
```
"""
function containment_problem(I::MPolyIdeal, n::Int; in_symbolic::Bool=false,
  use_minimal_primes::Bool=false,
  max_power::Int=typemax(Int))
  _check_ideal(I)
  @req n >= 1 "n must be positive"
  @req max_power >= 1 "max_power must be positive"
  @req !_is_unit_ideal(I) && !is_zero(I) "the ideal must be nonzero and proper"
  if !in_symbolic
    return _containment_threshold(I, n, use_minimal_primes, max_power,
      Dict{Int,typeof(I)}())
  end
  S = symbolic_power(I, n; use_minimal_primes)
  contained = if is_monomial(I)
    Se = _mono_exponents(S)
    k -> _mono_is_subset(Se, _mono_power_exponents(I, k))
  else
    k -> is_subset(S, I^k)
  end
  # I^(n) ⊆ I^k is monotone in k, true for k = 0, false for k > n.
  lo = use_minimal_primes ? 0 : div(n, big_height(I))
  hi = n
  while lo < hi
    mid = hi - div(hi - lo, 2) # Upper midpoint, without overflowing at the Int bound.
    if contained(mid)
      lo = mid
    else
      hi = mid - 1
    end
  end
  return lo
end
