@testset "intersection strategy, iterated Frobenius quotients, regularity" begin
  SP = SymbolicPowers

  # _intersect_all folds non-monomial lists and uses one backend call for
  # monomial lists. Compare with OSCAR's single intersection call.
  A, t = polynomial_ring(QQ, :t)
  N, _ = number_field(t^2 + 1, "i")
  for K in (QQ, GF(7), GF(2), N)
    R, (x, y, z, w) = polynomial_ring(K, [:x, :y, :z, :w])
    G, (a, b, c, d) = graded_polynomial_ring(K, [:a, :b, :c, :d];
      weights=[1, 2, 3, 1])
    lists = [
      [ideal(R, [x - y, z]), ideal(R, [x + z, w^2]), ideal(R, [y - w, x * z - 1]),
        ideal(R, [x^2 - y * z, w])],
      [ideal(R, [x, y])^2, ideal(R, [z, w])^2, ideal(R, [x, w])^3, ideal(R, [y, z])],
      [ideal(R, [x, y])^2, ideal(R, [x - z, y + w]), ideal(R, [x * y, z^3])],
      [ideal(R, [x + 1]), ideal(R, [one(R)]), ideal(R, [y^2 - 1, z])],
      [ideal(R, elem_type(R)[]), ideal(R, [x, y]), ideal(R, [z - w, x])],
      [ideal(R, [zero(R), x * y]), ideal(R, [y, z]), ideal(R, [x, w])],
      [ideal(R, [x - y]), ideal(R, [z^2 + w])],
      [ideal(G, [a^2 - b, c - a * b]), ideal(G, [d, b]), ideal(G, [a - d, c])],
      [ideal(G, [a * b, c]), ideal(G, [d^3, b^2]), ideal(G, [a, d])],
    ]
    for L in lists, k in 1:length(L)
      @test SP._intersect_all(L[1:k]) == (k == 1 ? L[1] : intersect(L[1:k]...))
    end
    @test_throws ArgumentError SP._intersect_all(typeof(ideal(R, [x]))[])
  end

  # Fermat cubic: the complete-intersection, primary and saturation methods
  # all end with intersections of many non-monomial ideals.
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  F3 = ideal(R, [x * (y^3 - z^3), y * (z^3 - x^3), z * (x^3 - y^3)])
  for n in 1:3
    S = symbolic_power(F3, n; method=:primary)
    @test symbolic_power(F3, n; method=:complete_intersection) == S
    @test symbolic_power(F3, n; method=:saturation) == S
    @test symbolic_power(F3, n; method=:primary, use_minimal_primes=true) == S
    pts = [P^n for P in minimal_primes(F3)]
    @test SP._intersect_all(pts) == intersect(pts...)
  end
  # embedded components selected in the primary method (several components)
  m = ideal(R, [x, y, z])
  E = intersect(ideal(R, [x * y, x * z, y * z]), m^3)
  for n in 1:3, mp in (false, true)
    comps = [
      Q for (Q, P) in primary_decomposition(E^n) if
      any(M -> is_subset(P, M), SP._defining_primes(E, mp))
    ]
    @test symbolic_power(E, n; method=:primary, use_minimal_primes=mp) ==
      intersect(comps...)
  end

  # Iterated quotients P^[q] : P : ... : P against the single colon by P^e.
  function single_colon(P, n)
    p = Int(characteristic(coefficient_ring(base_ring(P))))
    h = codim(P)
    q = p
    while q < n
      q *= p
    end
    e = h * (q - 1) - n + 1
    Pq = frobenius_power(P, q)
    return e <= 0 ? Pq : quotient(Pq, P^e)
  end
  for p in (2, 3, 5)
    F, (u, v, s) = polynomial_ring(GF(p), [:u, :v, :s])
    t345 = ideal(F, [v^2 - u * s, u^2 * v - s^2, u^3 - v * s])
    line = ideal(F, [u + v, s - 2 * u])
    point = ideal(F, [u - 1, v, s + u])
    curve = ideal(F, [u^2 - v * s + 1]) # principal: height 1
    for P in (t345, line, point, curve), n in 0:(p == 5 ? 5 : 4)
      expected = n == 0 ? ideal(F, [one(F)]) : (n == 1 ? P : single_colon(P, n))
      S = symbolic_power_frobenius(P, n)
      @test S == expected
      # The join needs no primary decomposition, which may not terminate
      # over these small prime fields.
      n >= 1 && @test S == symbolic_power_join(P, n; check=false)
    end
    @test symbolic_power_frobenius(ideal(F, elem_type(F)[]), 3) ==
      ideal(F, elem_type(F)[])
  end
  # graded ring, a nonprime finite field, and the twisted cubic
  K4, g = finite_field(2, 2, "g")
  G, (a, b, c, d) = graded_polynomial_ring(K4, [:a, :b, :c, :d])
  twc = ideal(G, [a * c - b^2, a * d - b * c, b * d - c^2])
  for n in 1:5
    @test symbolic_power_frobenius(twc, n) == single_colon(twc, n)
    @test symbolic_power(twc, n; method=:frobenius) ==
      symbolic_power_join(twc, n; check=false)
  end
  H, (a, b, c) = graded_polynomial_ring(GF(3), [:a, :b, :c]; weights=[3, 4, 5])
  P = ideal(H, [b^2 - a * c, a^2 * b - c^2, a^3 - b * c])
  for n in (2, 4, 5)
    @test symbolic_power_frobenius(P, n) == symbolic_power_join(P, n; check=false)
  end

  # Regularity samples from a minimal generating set agree with
  # cm_regularity on the generators returned by symbolic_power.
  function reference_samples(I, s)
    return [
      begin
        (S,) = SP._graded_copies(symbolic_power(I, n))
        r = Int(cm_regularity(S))
        (power=n, regularity=r, ratio=QQ(r, n))
      end for n in 1:s
    ]
  end
  pts = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1), (1, 2, 3)]
  ptideal(R, v, p) = ideal(R, [p[j] * v[i] - p[i] * v[j] for i in 1:3 for j in (i + 1):3])
  for K in (QQ, GF(7))
    R, v = polynomial_ring(K, [:x, :y, :z])
    I = ideal(R, gens(intersect([ptideal(R, v, p) for p in pts]...)))
    @test symbolic_regularity_samples(I; sample_size=2) == reference_samples(I, 2)
    G, (a, b, c, d) = graded_polynomial_ring(K, [:a, :b, :c, :d])
    twc = ideal(G, [a * c - b^2, a * d - b * c, b * d - c^2])
    @test symbolic_regularity_samples(twc; sample_size=3) == reference_samples(twc, 3)
    @test asymptotic_regularity(twc; sample_size=3, approximate=true) ==
      minimum(s.ratio for s in reference_samples(twc, 3))
    # nonhomogeneous given generators of a homogeneous ideal, unused variable
    S, (x, y, z, w) = polynomial_ring(K, [:x, :y, :z, :w])
    J = ideal(S, [x * y, x * z + x * y * w, y * z])
    @test symbolic_regularity_samples(J; sample_size=3) == reference_samples(J, 3)
    @test symbolic_regularity_samples(J; sample_size=2, method=:primary) ==
      reference_samples(J, 2)
  end
end
