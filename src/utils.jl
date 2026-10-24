###############################################################################
#  Internal helpers
###############################################################################

function _check_ideal(I::MPolyIdeal)
  R = base_ring(I)
  @req coefficient_ring(R) isa Field "the coefficient ring must be a field"
  return R
end

_check_power(n::Int) = @req n >= 0 "the exponent must be non-negative"

_unit_ideal(R::MPolyRing) = ideal(R, [one(R)])
_zero_ideal(R::MPolyRing) = ideal(R, elem_type(R)[])
_maximal_homogeneous_ideal(R::MPolyRing) = ideal(R, gens(R))

_is_unit_ideal(I::MPolyIdeal) = is_one(I)

# Remove duplicates from a list of ideals (ideal equality is mathematical,
# so `unique` cannot be used).
function _unique_ideals(L::Vector{T}) where {T<:MPolyIdeal}
  res = T[]
  for P in L
    any(Q -> Q == P, res) || push!(res, P)
  end
  return res
end

# A minimal primary decomposition has distinct radicals and no redundant
# components. OSCAR documents an exception over number fields: its backend can
# return redundant components. Their radicals need not be associated primes.
# Normalize only in that case, and cache the result without changing OSCAR's
# own decomposition. Intersections of primary ideals with the same radical
# remain primary.
function _irredundant_primary_decomposition(pd::Vector{Tuple{T,T}}) where {T<:MPolyIdeal}
  result = Tuple{T,T}[]
  for (Q, P) in pd
    _is_unit_ideal(Q) && continue
    i = findfirst(pair -> pair[2] == P, result)
    if isnothing(i)
      push!(result, (Q, P))
    else
      result[i] = (intersect(result[i][1], Q), P)
    end
  end
  i = 1
  while length(result) > 1 && i <= length(result)
    others = T[result[j][1] for j in eachindex(result) if j != i]
    if is_subset(_intersect_all(others), result[i][1])
      deleteat!(result, i)
    else
      i += 1
    end
  end
  return result
end

function _primary_decomposition(I::T) where {T<:MPolyIdeal}
  pd = primary_decomposition(I)
  coefficient_ring(base_ring(I)) isa Oscar.NumField || return pd
  return get_attribute!(
    I, :symbolic_powers_irredundant_primary_decomposition
  ) do
    _irredundant_primary_decomposition(pd)
  end::Vector{Tuple{T,T}}
end

# Distinct associated primes of I. Both the backend decomposition and its
# number-field normalization are cached on the ideal object.
function _associated_primes(I::MPolyIdeal)
  pd = _primary_decomposition(I)
  return _unique_ideals([P for (_, P) in pd])
end

# Elements of L that are maximal (resp. minimal) with respect to inclusion.
function _maximal_elements(L::Vector{T}) where {T<:MPolyIdeal}
  return T[P for P in L if !any(Q -> !(Q == P) && is_subset(P, Q), L)]
end

function _minimal_elements(L::Vector{T}) where {T<:MPolyIdeal}
  return T[P for P in L if !any(Q -> !(Q == P) && is_subset(Q, P), L)]
end

# The primes P used in the definition of the symbolic powers: the maximal
# associated primes (default) or the minimal primes. A primary component of a
# power of I survives in the symbolic power iff its prime lies in one of them.
function _defining_primes(I::MPolyIdeal, use_minimal_primes::Bool)
  primes = _associated_primes(I)
  return use_minimal_primes ? _minimal_elements(primes) : _maximal_elements(primes)
end

# Intersection of a nonempty list of ideals.
function _intersect_all(L::Vector{T}) where {T<:MPolyIdeal}
  @req !isempty(L) "cannot intersect an empty list of ideals"
  length(L) == 1 && return L[1]
  return intersect(L...)
end

###############################################################################
#  Degrees and homogeneity
###############################################################################

# Does the ring carry a positive Z-grading? Ungraded rings are regarded as
# standard graded.
function _has_positive_z_grading(R::MPolyRing)
  R isa MPolyDecRing || return true
  is_z_graded(R) || return false
  return all(x -> degree(Int, x) > 0, gens(R))
end

# Positive integer weights of the variables: the degrees in a Z-graded ring,
# all ones otherwise.
function _degree_weights(R::MPolyRing)
  if R isa MPolyDecRing
    @req is_z_graded(R) "degree invariants require a positive Z-grading"
    w = [Int(degree(Int, x)) for x in gens(R)]
    @req all(>(0), w) "all variables must have positive degree"
    return w
  end
  return ones(Int, ngens(R))
end

# Weight of an exponent vector, with overflow checks.
function _term_weight(e, w::Vector{Int})
  d = 0
  for i in eachindex(w)
    d = Base.Checked.checked_add(d, Base.Checked.checked_mul(w[i], e[i]))
  end
  return d
end

# Weighted degree of a polynomial (the maximum over its terms).
function _weighted_degree(f::MPolyRingElem, w::Vector{Int})
  @req !iszero(f) "the zero polynomial has no degree"
  return maximum(_term_weight(e, w) for e in exponents(f))
end

function _is_weighted_homogeneous(f::MPolyRingElem, w::Vector{Int})
  d = -1
  for e in exponents(f)
    t = _term_weight(e, w)
    if d < 0
      d = t
    elseif t != d
      return false
    end
  end
  return true
end

# Homogeneous generators of I, or `nothing` if I is not homogeneous. The given
# generators are tried first. Otherwise the reduced Gröbner basis decides: for
# a homogeneous ideal (positive weights) the reduced Gröbner basis with respect
# to any monomial ordering consists of homogeneous polynomials.
function _homogeneous_generators_or_nothing(I::MPolyIdeal)
  w = _degree_weights(base_ring(I))
  g = filter(!iszero, gens(I))
  all(f -> _is_weighted_homogeneous(f, w), g) && return g
  g = filter(!iszero, gens(groebner_basis(I; complete_reduction=true)))
  return all(f -> _is_weighted_homogeneous(f, w), g) ? g : nothing
end

function _homogeneous_generators(I::MPolyIdeal)
  g = _homogeneous_generators_or_nothing(I)
  @req g !== nothing "the ideal must be homogeneous"
  return g
end

_is_homogeneous_ideal(I::MPolyIdeal) = _homogeneous_generators_or_nothing(I) !== nothing

# A graded polynomial ring T on the variables of R (R itself if it is already
# graded, otherwise the standard graded ring) and the identification R -> T.
function _graded_ring(R::MPolyRing)
  R isa MPolyDecRing && return R, identity
  T, _ = graded_polynomial_ring(coefficient_ring(R), symbols(R); cached=false)
  return T, hom(R, T, gens(T))
end

# The homogeneous ideals in `Is` (all in the same ring), moved to a common
# graded ring.
function _graded_copies(Is::MPolyIdeal...)
  R = base_ring(first(Is))
  _degree_weights(R)
  T, to_T = _graded_ring(R)
  return map(Is) do I
    @req base_ring(I) === R "the ideals must live in the same ring"
    g = [to_T(f) for f in _homogeneous_generators(I)]
    isempty(g) ? ideal(T, elem_type(T)[]) : ideal(T, g)
  end
end

###############################################################################
#  Monomial ideals
###############################################################################

_is_monomial_poly(f::MPolyRingElem) = length(f) == 1

# Minimal monomial generators (monic) of a monomial ideal. If the given
# generators are monomials no Gröbner basis is needed.
function _minimal_monomial_generators(I::MPolyIdeal)
  R = base_ring(I)
  @req is_monomial(I) "the ideal must be a monomial ideal"
  is_zero(I) && return elem_type(R)[]
  supplied = filter(!iszero, gens(I))
  G = if all(_is_monomial_poly, supplied)
    supplied
  else
    gens(groebner_basis(I; complete_reduction=true))
  end
  mons = unique([monomial(f, 1) for f in G if !iszero(f)])
  exps = [exponent_vector(m, 1) for m in mons]
  keep = [
    i for i in eachindex(mons) if
    !any(j -> exps[j] != exps[i] && all(exps[j] .<= exps[i]), eachindex(mons))
  ]
  return sort(mons[keep]; by=m -> (total_degree(m), exponent_vector(m, 1)))
end

function _is_squarefree_monomial_ideal(I::MPolyIdeal)
  is_monomial(I) || return false
  return all(m -> all(<=(1), exponent_vector(m, 1)), _minimal_monomial_generators(I))
end

# Supports of the minimal monomial generators, as sets of variable indices.
_monomial_supports(I::MPolyIdeal) =
  [BitSet(findall(>(0), exponent_vector(m, 1))) for m in _minimal_monomial_generators(I)]

###############################################################################
#  Frobenius exponents
###############################################################################

function _check_frobenius_exponent(R::MPolyRing, q::Int)
  @req q >= 1 "the exponent must be positive"
  q == 1 && return nothing
  p = characteristic(coefficient_ring(R))
  @req p > 0 "a nontrivial Frobenius power requires positive characteristic"
  a = ZZ(q)
  while a > 1 && is_zero(mod(a, p))
    a = div(a, p)
  end
  @req a == 1 "the exponent must be a power of the characteristic"
  return nothing
end
