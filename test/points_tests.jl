# Dimension-one homogeneous ideals (points in weighted projective space): the
# intersection of powers of rational point ideals, the complete-intersection
# shortcut and the fallbacks are compared with the general primary method or
# with the previous saturation code path.

# The prime of the point p (homogeneous coordinates) in the polynomial ring R.
function _point_ideal(R, p)
  x = gens(R)
  j = findfirst(!iszero, p)
  return ideal(R, [p[j] * x[i] - p[i] * x[j] for i in eachindex(x) if i != j])
end

_points_ideal(R, pts) = ideal(R, gens(intersect([_point_ideal(R, p) for p in pts]...)))

# Outcome of the search for the points by projections (:ok, :fail or
# :nonrational), on a fresh copy of I.
function _projection_status(I)
  J = ideal(base_ring(I), gens(I))
  pairs = SymbolicPowers._equal_weight_pairs(base_ring(J))
  return SymbolicPowers._project_points!(typeof(J)[], J, pairs, 1)
end

# The code path of method :saturation before the point shortcut: I^n if the
# irrelevant ideal m is associated to I (default convention), otherwise the
# saturation of I^n with respect to m.
function _saturation_reference(R, gs, n, mp)
  n == 0 && return ideal(R, [one(R)])
  I = ideal(R, gs)
  m = ideal(R, gens(R))
  !mp && !(saturation(I, m) == I) && return I^n
  return saturation(I^n, m)
end

# Compare :saturation and :auto with the general method (`reference = :primary`)
# or with the saturation path above on fresh ideal objects, for both
# conventions and several exponents (a repeated exponent reuses the caches).
# For radical ideals both conventions agree, so one reference suffices.
function _check_points(R, gs; ns=(0, 1, 2, 3, 2), radical=false, reference=:primary)
  I = ideal(R, gs)
  refs = Dict{Tuple{Int,Bool},typeof(I)}()
  for n in ns, mp in (false, true)
    ref = get!(refs, (n, radical ? false : mp)) do
      if reference == :primary
        symbolic_power(ideal(R, gs), n; method=:primary, use_minimal_primes=mp)
      else
        _saturation_reference(R, gs, n, mp)
      end
    end
    @test symbolic_power(I, n; method=:saturation, use_minimal_primes=mp) == ref
    @test symbolic_power(ideal(R, gs), n; use_minimal_primes=mp) == ref
  end
  return I
end

@testset "points: rational point ideals" begin
  for K in (QQ, GF(2), GF(32003))
    R, (x, y, z) = polynomial_ring(K, [:x, :y, :z])
    pts = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1], [1, 1, 0], [0, 1, 1]]
    K == QQ && append!(pts, [[1, 2, 3], [3, 1, -2]])
    I = _check_points(R, gens(_points_ideal(R, pts)); radical=true)
    # The shortcut was taken, and its data are cached on the ideal.
    @test SymbolicPowers._rational_point_components(I) isa Vector
    # found by projections, also where (1:1:1) and (0:1:1) have the same image
    @test _projection_status(I) == :ok
    @test !SymbolicPowers._irrelevant_ideal_is_associated(I)
    @test symbolic_power(I, 4) == intersect([_point_ideal(R, p)^4 for p in pts]...)
  end
  # Nemo's native prime fields
  F7 = Nemo.Native.GF(7)
  R7, _ = polynomial_ring(F7, [:x, :y, :z])
  I = _check_points(R7,
    gens(_points_ideal(R7, [[1, 0, 0], [0, 1, 0], [1, 1, 1], [1, 2, 3], [3, 1, 5]]));
    ns=(2, 3), radical=true)
  @test SymbolicPowers._rational_point_components(I) isa Vector
  # number field coordinates: the saturation is used (a performance choice)
  Qt, t = polynomial_ring(QQ, :t)
  L, r = number_field(t^2 - 2, :r)
  T, (x, y, z) = polynomial_ring(L, [:x, :y, :z])
  I = _check_points(T,
    gens(_points_ideal(T, [[1, r, 0], [1, -r, 0], [0, 1, 1], [1, 1, r]]));
    ns=(2, 3), radical=true, reference=:saturation)
  @test SymbolicPowers._rational_point_components(I) === nothing
  # GF(9): the points are found by projections (successive hyperplane
  # sections also reduce fat points). If the projections give no answer (for
  # non-standard weights), the saturation is used.
  F9, a = finite_field(3, 2, "a")
  S, _ = polynomial_ring(F9, [:x, :y, :z])
  F9pts = [[1, a, 0], [1, a + 1, a], [0, 1, a^2], [1, 1, 1]]
  I = _check_points(S, gens(_points_ideal(S, F9pts)); ns=(2, 3), radical=true)
  @test SymbolicPowers._rational_point_components(I) isa Vector
  F9fat = intersect(_point_ideal(S, F9pts[1])^2, _point_ideal(S, F9pts[2]))
  I = _check_points(S, gens(F9fat); ns=(2,))
  @test SymbolicPowers._rational_point_components(I) isa Vector
  W9, (u, v, w) = graded_polynomial_ring(F9, [:u, :v, :w]; weights=[1, 2, 3])
  I = intersect(ideal(W9, [v, w]), ideal(W9, [u, w]), ideal(W9, [v - a * u^2, w - u^3]))
  @test _projection_status(I) == :fail
  I = _check_points(W9, gens(I); ns=(2,), radical=true)
  @test SymbolicPowers._rational_point_components(I) === nothing
  # points in P^3, collinear points (a complete intersection), a grid
  R4, _ = polynomial_ring(QQ, [:a, :b, :c, :d])
  I = _check_points(R4,
    gens(
      _points_ideal(
        R4, [[1, 0, 0, 0], [0, 1, 0, 0], [1, 1, 1, 1],
          [1, 2, 3, 4], [2, -1, 1, 3]]
      ),
    ); ns=(2, 3), radical=true)
  @test !SymbolicPowers._is_points_complete_intersection(ideal(R4, gens(I)))
  # complete intersections: I^(n) = I^n without further computation
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  I = _check_points(R, gens(_points_ideal(R, [[1, i, 0] for i in 0:4])); radical=true)
  @test SymbolicPowers._is_points_complete_intersection(I)
  grid = _points_ideal(R, [[i, j, 1] for i in 0:1 for j in 0:2])
  @test length(gens(grid)) == 2
  I = _check_points(R, gens(grid); ns=(2, 3), radical=true)
  @test SymbolicPowers._is_points_complete_intersection(I)
  @test symbolic_power(I, 3; method=:saturation) == grid^3
  # only the given generators are counted: with a redundant generator the
  # points are found by projections
  I = _check_points(R, [gens(grid); sum(gens(grid))]; ns=(2, 3), radical=true)
  @test !SymbolicPowers._is_points_complete_intersection(I)
  @test SymbolicPowers._rational_point_components(I) isa Vector
  @test symbolic_power(I, 3; method=:saturation) == grid^3
  # a fat complete intersection (non-radical), in a standard graded ring
  G, (u, v, w) = graded_polynomial_ring(GF(7), [:u, :v, :w])
  I = _check_points(G, [u^2 * (u - w), v^3 - u * w^2]; ns=(2, 3))
  @test SymbolicPowers._is_points_complete_intersection(I)
end

@testset "points: gradings and generators" begin
  # standard graded ring
  G, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
  I = _points_ideal(G, [[1, 0, 0], [0, 1, 1], [1, 2, 3], [1, 1, 1]])
  _check_points(G, gens(I); ns=(1, 2, 3), radical=true)
  # inhomogeneous generators of a homogeneous ideal
  R, (u, v, w) = polynomial_ring(QQ, [:u, :v, :w])
  g = gens(_points_ideal(R, [[1, 0, 0], [0, 1, 0], [1, 1, 1], [1, 2, 4]]))
  inhom = [g[1]; [g[i] + u * g[1] for i in 2:length(g)]]
  @test !all(f -> SymbolicPowers._is_weighted_homogeneous(f, [1, 1, 1]), inhom)
  _check_points(R, inhom; ns=(2, 3), radical=true)
  # weighted P(1,1,2): points whose ideals are generated by linear forms, and
  # points whose ideals are not (s - u^2 is not linear), which fall back
  W, (u, v, s) = graded_polynomial_ring(QQ, [:u, :v, :s]; weights=[1, 1, 2])
  lin = intersect(ideal(W, [u - v, s]), ideal(W, [u, v]), ideal(W, [u + 2 * v, s]),
    ideal(W, [u - 3 * v, s]))
  I = _check_points(W, gens(lin); radical=true)
  @test SymbolicPowers._rational_point_components(I) isa Vector
  nonlin = intersect(ideal(W, [u - v, s - u^2]), ideal(W, [u, v]),
    ideal(W, [u + 2 * v, s - 3 * v^2]), ideal(W, [v, s - u^2]))
  I = _check_points(W, gens(nonlin); ns=(2, 3), radical=true)
  @test SymbolicPowers._rational_point_components(I) === nothing
end

@testset "points: fallbacks and embedded components" begin
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  P = [_point_ideal(R, p) for p in ([1, 0, 0], [0, 1, 0], [1, 1, 1], [1, 2, 3])]
  # fat points: I^(n) is the intersection of the ordinary powers P_i^{n k_i}
  fat = intersect(P[1]^2, P[2], P[3], P[4])
  I = _check_points(R, gens(fat); ns=(2, 3))
  # successive hyperplane sections reduce the double point
  @test _projection_status(I) == :ok
  comps = SymbolicPowers._rational_point_components(I)
  @test length(comps) == 4 && all(((Q, k),) -> k == (Q == P[1] ? 2 : 1), comps)
  fat2 = intersect(P[1]^3, P[2]^2, P[3], P[4]^2)
  I = _check_points(R, gens(fat2); ns=(2, 3))
  @test SymbolicPowers._rational_point_components(I) isa Vector
  @test symbolic_power(I, 2) == intersect(P[1]^6, P[2]^4, P[3]^2, P[4]^4)
  # a curvilinear double point (x^2, y) is not a power of its prime: fallback
  curv = intersect(ideal(R, [x^2, y]), P[2], P[3], P[4])
  I = _check_points(R, gens(curv); ns=(2, 3))
  @test SymbolicPowers._rational_point_components(I) === nothing
  # a Galois orbit (three conjugate points, not a complete intersection)
  orbit = ideal(R, [x^2 - y * z, x * y - 2 * z^2, y^2 - 2 * x * z])
  I = _check_points(R, gens(intersect(orbit, P[3])); ns=(2, 3), radical=true)
  @test _projection_status(I) == :nonrational
  @test SymbolicPowers._rational_point_components(I) === nothing
  # a conjugate pair with a rational projection to (y : z), detected by the
  # next projection
  I = _check_points(R, gens(intersect(ideal(R, [y - z, x^2 + z^2]), P[1], P[4]));
    ns=(2, 3), radical=true)
  @test _projection_status(I) == :nonrational
  @test SymbolicPowers._rational_point_components(I) === nothing
  # rational points with an embedded irrelevant component: the default
  # convention gives I^n, the minimal-prime convention the intersection of
  # the powers of the point ideals (the saturation of I is radical)
  m = ideal(R, gens(R))
  I = _check_points(R, gens(intersect(P..., m^3)); ns=(1, 2, 3))
  @test SymbolicPowers._irrelevant_ideal_is_associated(I)
  @test SymbolicPowers._rational_point_components(I) isa Vector
  @test symbolic_power(I, 3; method=:saturation) == I^3
  @test symbolic_power(I, 3; use_minimal_primes=true) == intersect([Q^3 for Q in P]...)
  # fat points with an embedded component
  I = _check_points(R, gens(intersect(fat, m^4)); ns=(2,))
  @test SymbolicPowers._irrelevant_ideal_is_associated(I)
  @test SymbolicPowers._rational_point_components(I) isa Vector
  # Candidates are certified: m, a prime not containing I, or a non-linear
  # prime are rejected.
  I = ideal(R, gens(intersect(P...)))
  @test SymbolicPowers._is_rational_point_prime(P[1], I)
  @test !SymbolicPowers._is_rational_point_prime(m, I)
  @test !SymbolicPowers._is_rational_point_prime(_point_ideal(R, [1, 3, 3]), I)
  @test !SymbolicPowers._is_rational_point_prime(orbit, ideal(R, gens(orbit)))
  @test symbolic_power(I, 2; method=:saturation) == intersect([Q^2 for Q in P]...)
  # :auto uses the saturation method in dimension one; for the minimal-prime
  # convention only when a shortcut applies (rational points, complete
  # intersections), since saturating I^n can be much slower than the
  # componentwise method
  @test SymbolicPowers._choose_method(ideal(R, gens(fat)), true) == :saturation
  @test SymbolicPowers._choose_method(ideal(R, [x * y, x * z]), true) == :monomial
  @test SymbolicPowers._choose_method(ideal(R, [x * y - z^2, x - y]), true) == :saturation
  pair = intersect(ideal(R, [y - z, x^2 + z^2]), P[1], P[4])
  @test SymbolicPowers._choose_method(ideal(R, gens(pair)), true) == :components
  @test SymbolicPowers._choose_method(ideal(R, gens(pair)), false) == :saturation
  @test symbolic_power(ideal(R, gens(pair)), 2; use_minimal_primes=true) ==
    symbolic_power(ideal(R, gens(pair)), 2; method=:primary, use_minimal_primes=true)
  S, (a, b, c, d) = polynomial_ring(QQ, [:a, :b, :c, :d])
  @test SymbolicPowers._choose_method(ideal(S, [a * b - c * d]), true) == :unmixed
end

@testset "points: projective line and other weights" begin
  # points of P^1 (principal ideals, I^(n) = I^n), rational or not
  for K in (QQ, GF(3))
    S, (u, v) = polynomial_ring(K, [:u, :v])
    _check_points(S, [u * v * (u - v)]; ns=(2, 3), radical=true)
    _check_points(S, [u^2 * v * (u + v)^3]; ns=(2, 3))
    _check_points(S, [u * (u^2 + v^2)]; ns=(2,), radical=true)
  end
  # P(1,2,3): no two variables of equal weight, so the projections give no
  # answer and the saturation is used
  W, (a, b, c) = graded_polynomial_ring(QQ, [:a, :b, :c]; weights=[1, 2, 3])
  pts = intersect(ideal(W, [b, c]), ideal(W, [a, b]), ideal(W, [a, c]),
    ideal(W, [b - a^2, c - a^3]))
  @test _projection_status(pts) == :fail
  I = _check_points(W, gens(pts); ns=(2, 3), radical=true)
  @test SymbolicPowers._rational_point_components(I) === nothing
  # a single rational point is recognized without projections
  @test SymbolicPowers._rational_point_components(ideal(W, [a, b])) ==
    [(ideal(W, [a, b]), 1)]
  # a weighted fat point that is not a complete intersection
  _check_points(W, [a^4, a^2 * b, b^2]; ns=(2, 3))
  # unused variable: the points lie in the hyperplane d = 0 of P^3
  R4, (x, y, z, d) = polynomial_ring(GF(5), [:x, :y, :z, :d])
  I = _check_points(R4,
    gens(
      _points_ideal(
        R4, [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [1, 1, 1, 0], [1, 2, 3, 0]]
      ),
    );
    ns=(2, 3), radical=true)
  @test SymbolicPowers._rational_point_components(I) isa Vector
end
