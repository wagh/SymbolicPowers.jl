###############################################################################
#  Unmixedness without primary decomposition; componentwise symbolic powers
###############################################################################

# Are the associated primes of I already available without a new primary
# decomposition (package cache or OSCAR's cached decomposition)?
_associated_primes_known(I::MPolyIdeal) =
  has_attribute(I, :symbolic_powers_associated_primes) ||
  Oscar.is_known(primary_decomposition, I)

# Whether all associated primes of I have the same height, i.e. whether
# bight(I) == codim(I). The equidimensional hull of I is the intersection of the
# primary components of I of maximal dimension (OSCAR, Singular's equidimMax),
# so I ⊆ hull(I) always, with equality iff every associated prime has
# dimension dim(R/I), i.e. height codim(I). This avoids a primary decomposition
# of I; if the associated primes are already known they are used instead.
# The value is exponent independent and cached on the ideal.
function _is_unmixed(I::MPolyIdeal)
  return get_attribute!(I, :symbolic_powers_is_unmixed) do
    if _associated_primes_known(I)
      c = codim(I)
      return all(P -> codim(P) == c, _associated_primes(I))
    end
    Oscar.is_known(is_prime, I) && is_prime(I) && return true
    return is_subset(equidimensional_hull(I), I)
  end::Bool
end

# Cheap unmixedness shortcut for invariants of nonmonomial ideals: if it
# returns true, every associated prime has height codim(I). (Monomial ideals
# keep their own combinatorial paths.)
_unmixed_shortcut(I::MPolyIdeal) = !is_monomial(I) && _is_unmixed(I)

# Componentwise method, used by `:auto` for nonmonomial ideals that are not
# unmixed. Let I = ∩ Q_P be the (irredundant) primary decomposition of I and
# let M be a defining prime (a maximal associated prime, or a minimal prime for
# the minimal-prime convention). Localizing at M kills the components whose
# prime is not contained in M, so
#
#   I^n R_M ∩ R = J^n R_M ∩ R   for every ideal J with J R_M = I R_M.       (*)
#
# 1. Defining primes M that are minimal primes of I ("isolated" targets): only
#    Q_M survives in R_M. Group these primes by height h and put
#    I_h = ∩_{ht M = h} Q_M. The minimal primes of I_h are exactly these M (all
#    of height h, pairwise incomparable), and I_h R_M = Q_M R_M = I R_M. Every
#    associated prime of I_h^n contains one of them, so the primary components
#    of I_h^n of maximal dimension are exactly the components at these M, i.e.
#    equidimensional_hull(I_h^n) = ∩_{ht M = h} (I^n R_M ∩ R) by (*).
# 2. The remaining (embedded) defining primes: put J = ∩ Q_P over the primes P
#    contained in one of them. For such an M, J R_M = I R_M, so (*) and the
#    selection rule of the primary method applied to a primary decomposition of
#    J^n give their contribution.
# Then I^(n) is the intersection of these parts. If no defining prime is
# isolated, J = I and this is the primary method. Otherwise the potentially
# hard primary decomposition of I^n is replaced by hulls of powers of the
# isolated parts (and a decomposition of the smaller J^n). For the minimal-prime
# convention, and for ideals without embedded primes (e.g. radical ideals of
# mixed height), every defining prime is isolated and no power is decomposed.
# The argument holds over every coefficient field; it only uses the
# decomposition of I, as the primary method does.
#
# Not used: for an unmixed ideal of height h and an element f of the ideal of
# h x h minors of the Jacobian with codim(I + (f)) > h, one has I^(n) =
# I^n : f^∞ (the Jacobian criterion makes I_Q a complete intersection in a
# regular system of parameters at every prime Q not containing f). In OSCAR
# this saturation was at most about 4x faster than equidimensional_hull(I^n) on
# small primes, but much slower whenever f is not sparse (three skew lines in
# P^3, n = 6: 244 s against 0.11 s), so the hulls are kept.
function _symbolic_power_components(
  I::T, n::Int, use_minimal_primes::Bool
) where {T<:MPolyIdeal}
  n == 1 && return use_minimal_primes ? minimal_part(I) : I
  minimal = _defining_primes(I, true)
  targets = _defining_primes(I, use_minimal_primes)
  isolated = T[M for M in targets if any(P -> P == M, minimal)]
  isempty(isolated) && return _symbolic_power_primary(I, n, use_minimal_primes)
  embedded = T[M for M in targets if !any(P -> P == M, isolated)]
  pd = _primary_decomposition(I)
  parts = T[]
  heights = [codim(M) for M in isolated]
  for h in unique(heights)
    primes_h = isolated[heights .== h]
    I_h = _intersect_all(T[Q for (Q, P) in pd if any(M -> M == P, primes_h)])
    push!(parts, equidimensional_hull(I_h^n))
  end
  if !isempty(embedded)
    J = _intersect_all(T[Q for (Q, P) in pd if any(M -> is_subset(P, M), embedded)])
    push!(parts, _symbolic_power_primary_targets(J^n, embedded))
  end
  return _intersect_all(parts)
end

# The intersection of the primary components of `power` whose prime lies in
# one of `targets` (see `_symbolic_power_primary`).
function _symbolic_power_primary_targets(power::T, targets::Vector{T}) where {T<:MPolyIdeal}
  comps = T[
    Q for (Q, P) in primary_decomposition(power) if any(M -> is_subset(P, M), targets)
  ]
  return _intersect_all(comps)
end
