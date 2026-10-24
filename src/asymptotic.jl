###############################################################################
#  Symbolic defect
###############################################################################

@doc raw"""
    symbolic_defect(I::MPolyIdeal, n::Int; kwargs...) -> Int

Return the symbolic defect of `I` in degree `n`, the minimal number of
generators of ``I^{(n)} / I^n``. The ideal must be homogeneous. Keyword
arguments are passed to [`symbolic_power`](@ref).

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> symbolic_defect(ideal(R, [x*y, x*z, y*z]), 2)
1
```
"""
function symbolic_defect(I::MPolyIdeal, n::Int; kwargs...)
  _check_ideal(I)
  _check_power(n)
  (n == 0 || is_zero(I) || _is_unit_ideal(I)) && return 0
  S = symbolic_power(I, n; kwargs...)
  Ig, Sg = _graded_copies(I, S)
  T = base_ring(Ig)
  A, pr = quo(T, Ig^n)
  J = ideal(A, [pr(f) for f in gens(Sg)])
  return length(minimal_generating_set(J))
end

###############################################################################
#  Symbolic polyhedron and Waldschmidt constant
###############################################################################

# Newton polyhedron conv(exponents) + R_{\ge 0}^d of a monomial ideal.
function _newton_polyhedron(Q::MPolyIdeal)
  d = ngens(base_ring(Q))
  E = [exponent_vector(m, 1) for m in _minimal_monomial_generators(Q)]
  V = matrix(QQ, reduce(vcat, [permutedims(e) for e in E]))
  rays = identity_matrix(QQ, d)
  return convex_hull(V, rays)
end

@doc raw"""
    symbolic_polyhedron(I::MPolyIdeal;
                        use_minimal_primes::Bool = false) -> Polyhedron

Return the symbolic polyhedron of the monomial ideal `I`: the intersection of
the Newton polyhedra of the ideals ``Q_{\subseteq P}``, where ``P`` runs over
the maximal associated primes of `I` (Cooper–Embree–Hà–Hoefel).
Exponents of monomials in ``I^{(n)}`` lie in the dilation ``nP``. The converse
can fail: `(x^2,y^2)` has `(1,1)` in its Newton polyhedron but `x*y` is not
in the ideal. For squarefree monomial ideals the converse does hold.
With `use_minimal_primes = true`, only minimal primary components are used.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> P = symbolic_polyhedron(ideal(R, [x*y, y*z, x*z]));

julia> sort(Vector{QQFieldElem}.(collect(vertices(P)))) ==
       sort([QQ.([0, 1, 1]), QQ.([1, 0, 1]), QQ.([1, 1, 0]),
             [QQ(1, 2), QQ(1, 2), QQ(1, 2)]])
true
```
"""
function symbolic_polyhedron(I::MPolyIdeal; use_minimal_primes::Bool=false)
  _check_ideal(I)
  @req is_monomial(I) "the symbolic polyhedron is only defined for monomial ideals"
  @req !_is_unit_ideal(I) && !is_zero(I) "the ideal must be nonzero and proper"
  polys = [
    _newton_polyhedron(Q) for Q in _localized_components(I, use_minimal_primes)
  ]
  length(polys) == 1 && return polys[1]
  return intersect(polys...)
end

@doc raw"""
    waldschmidt_constant(I::MPolyIdeal;
                        use_minimal_primes::Bool = false) -> QQFieldElem

Return the Waldschmidt constant
``\widehat\alpha(I) = \lim_{n\to\infty} \alpha(I^{(n)})/n``.

The ideal must be nonzero, proper, and monomial. The answer is exact: minimize
the positive (weighted) degree on the symbolic polyhedron. For other
homogeneous ideals use [`waldschmidt_upper_bound`](@ref); finite sampling
does not certify an exact constant.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> waldschmidt_constant(ideal(R, [x*y, y*z, x*z]))
3//2
```
"""
function waldschmidt_constant(I::MPolyIdeal; use_minimal_primes::Bool=false)
  R = _check_ideal(I)
  @req !_is_unit_ideal(I) && !is_zero(I) "the ideal must be nonzero and proper"
  @req is_monomial(I) "exact Waldschmidt constants are implemented for monomial ideals; use waldschmidt_upper_bound for samples"
  w = _degree_weights(R)
  N = symbolic_polyhedron(I; use_minimal_primes)
  return minimum(sum(QQ(w[i]) * v[i] for i in eachindex(w)) for v in vertices(N))
end

@doc raw"""
    waldschmidt_upper_bound(I::MPolyIdeal; sample_size::Int = 5,
                           kwargs...) -> QQFieldElem

Return the certified upper bound
``\min_{1 \le n \le \mathrm{sample\_size}} \alpha(I^{(n)})/n`` for a nonzero
proper homogeneous ideal in a positive Z-grading (the standard grading for
an undecorated ring). The symbolic powers form a multiplicative filtration,
so the sequence of initial degrees is subadditive and Fekete's lemma applies.
An equality with the limit is not asserted. Keywords in `kwargs` are passed
to [`symbolic_power`](@ref).

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> waldschmidt_upper_bound(ideal(R, [x*y, x*z, y*z]); sample_size=1)
2
```
"""
function waldschmidt_upper_bound(I::MPolyIdeal; sample_size::Int=5, kwargs...)
  _check_ideal(I)
  @req !_is_unit_ideal(I) && !is_zero(I) "the ideal must be nonzero and proper"
  @req sample_size >= 1 "sample_size must be positive"
  _homogeneous_generators(I)
  return minimum(
    QQ(initial_degree(symbolic_power(I, n; kwargs...)), n) for n in 1:sample_size
  )
end

###############################################################################
#  Resurgence and asymptotic regularity
###############################################################################

@doc raw"""
    lower_bound_resurgence(I::MPolyIdeal; sample_size::Int = 5,
                           use_waldschmidt::Bool = false) -> QQFieldElem

Return a lower bound for the resurgence
``\rho(I) = \sup\{m/n : I^{(m)} \not\subseteq I^n\}``, namely the maximum of 1
and ``(m_n - 1)/n`` for ``2 \le n \le \text{sample\_size}``, where ``m_n`` is
the smallest `m` with ``I^{(m)} \subseteq I^n``. With
`use_waldschmidt = true`, also use ``\alpha(I)/U``, where `U` is the exact
monomial Waldschmidt constant or a certified sampled upper bound for other
homogeneous ideals. Dividing by an upper bound still gives a lower bound.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> lower_bound_resurgence(ideal(R, [x*y, x*z, y*z]))
6//5
```
"""
function lower_bound_resurgence(I::MPolyIdeal; sample_size::Int=5,
  use_waldschmidt::Bool=false)
  _check_ideal(I)
  @req !_is_unit_ideal(I) && !is_zero(I) "the ideal must be nonzero and proper"
  @req sample_size >= 1 "sample_size must be positive"
  # symbolic powers are shared between the containment searches
  cache = Dict{Int,typeof(I)}()
  bound = QQ(1)
  for n in 2:sample_size
    m = _containment_threshold(I, n, false, typemax(Int), cache)
    bound = max(bound, QQ(m - 1, n))
  end
  if use_waldschmidt
    upper = if is_monomial(I)
      waldschmidt_constant(I)
    else
      _homogeneous_generators(I)
      minimum(
        QQ(initial_degree(_cached_symbolic_power!(cache, I, n, false)), n)
        for n in 1:sample_size
      )
    end
    bound = max(bound, QQ(initial_degree(I)) / upper)
  end
  return bound
end

@doc raw"""
    symbolic_regularity_samples(I::MPolyIdeal; sample_size::Int = 10,
                               kwargs...) -> Vector

Return named tuples `(power, regularity, ratio)` for symbolic powers 1 through
`sample_size`. Only nonzero proper ideals homogeneous in the standard grading
are accepted. Ratios are exact rationals. Finite samples neither determine an
asymptotic limit nor provide a general upper or lower bound for it.
Keywords in `kwargs` are passed to [`symbolic_power`](@ref).
"""
function symbolic_regularity_samples(I::MPolyIdeal; sample_size::Int=10, kwargs...)
  R = _check_ideal(I)
  @req !_is_unit_ideal(I) && !is_zero(I) "the ideal must be nonzero and proper"
  @req sample_size >= 1 "sample_size must be positive"
  @req all(==(1), _degree_weights(R)) "regularity requires the standard Z-grading"
  _homogeneous_generators(I)
  vals = NamedTuple{(:power, :regularity, :ratio),Tuple{Int,Int,QQFieldElem}}[]
  for n in 1:sample_size
    (S,) = _graded_copies(symbolic_power(I, n; kwargs...))
    r = Int(cm_regularity(S))
    push!(vals, (power=n, regularity=r, ratio=QQ(r, n)))
  end
  return vals
end

@doc raw"""
    asymptotic_regularity(I::MPolyIdeal; sample_size::Int = 10,
                         approximate::Bool = false, kwargs...) -> QQFieldElem

Legacy finite-sampling interface. Requires `approximate = true` to explicitly
request the minimum of the sampled ratios in [`symbolic_regularity_samples`](@ref).
The returned number is a finite statistic, not a certified asymptotic
regularity or a bound on it. Prefer `symbolic_regularity_samples`, which also
reports every sampled power and regularity. Without explicit opt-in an
`ArgumentError` is thrown.
"""
function asymptotic_regularity(I::MPolyIdeal; sample_size::Int=10,
  approximate::Bool=false, kwargs...)
  @req approximate "finite sampling does not compute asymptotic regularity; use symbolic_regularity_samples or pass approximate=true"
  return minimum(s.ratio for s in symbolic_regularity_samples(I; sample_size, kwargs...))
end
