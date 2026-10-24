@testset "number-field decompositions and redundant backend components" begin
  A, t = polynomial_ring(QQ, :t)
  K, a = number_field(t^2 + 1, "a")
  R, (x, y, z) = polynomial_ring(K, [:x, :y, :z])
  I = ideal(R, [x * y, x * z, y * z])
  m = ideal(R, [x, y, z])
  raw = primary_decomposition(I)
  # OSCAR explicitly permits redundant components over number fields. Inject
  # a mathematically valid decomposition of that form, rather than depend on
  # whether a particular backend version happens to generate a redundancy.
  redundant = vcat(raw, [(m^2, m)])
  @test intersect(first.(redundant)...) == I
  set_attribute!(I, :primary_decomposition, redundant)
  @test big_height(I) == 2
  @test associated_prime_heights(I) == [2]
  @test !is_symbolic_equal_ordinary(I, 2)
  @test waldschmidt_constant(I) == QQ(3, 2)
  for mp in (false, true), n in 1:3
    expected = intersect(ideal(R, [x, y])^n, ideal(R, [x, z])^n, ideal(R, [y, z])^n)
    @test symbolic_power(I, n; method=:primary, use_minimal_primes=mp) == expected
    @test symbolic_power(I, n; method=:monomial, use_minimal_primes=mp) == expected
  end
  @test primary_decomposition(I) === redundant # leave OSCAR's cache untouched
  @test SymbolicPowers._primary_decomposition(I) ===
    SymbolicPowers._primary_decomposition(I)

  # Merge equal-radical primary components before taking their powers.
  P = ideal(R, [x, y])
  set_attribute!(P, :primary_decomposition, [(P, P), (P, P)])
  @test symbolic_power(P, 2; method=:complete_intersection) == P^2
  U = ideal(R, [one(R)])
  norm = SymbolicPowers._irredundant_primary_decomposition([(P, P), (P, P), (U, U)])
  @test length(norm) == 1 && norm[1] == (P, P)

  # A genuine embedded prime must survive normalization.
  H = ideal(R, [x^2, x * y])
  hpd = primary_decomposition(H)
  set_attribute!(H, :primary_decomposition, vcat(hpd, [(m^2, m)]))
  @test associated_prime_heights(H) == [1, 2]
  @test minimal_part(H) == ideal(R, [x])
  @test symbolic_power(H, 2; method=:primary) == H^2

  # Actual backend computations, with nonrational coefficients and without
  # injected decompositions, exercise the public number-field paths as well.
  Q = ideal(R, [x - a * y, z])
  for method in (:primary, :unmixed, :complete_intersection, :join)
    @test symbolic_power(Q, 2; method) == Q^2
  end
  G, (u, v, w) = graded_polynomial_ring(K, [:u, :v, :w])
  J = ideal(G, [u * v, u * w, v * w])
  @test symbolic_defect(J, 2) == 1
  @test symbolic_power(J, 2; method=:saturation) == J^2 + ideal(G, [u * v * w])
end

@testset "new associated primes contained in an old embedded prime" begin
  R, (x, y, z, w) = polynomial_ring(QQ, [:x, :y, :z, :w])
  J = ideal(R, [(x + y) * y, (x + y) * z, y * z])
  m = ideal(R, gens(R))
  I = intersect(J, m^3)
  @test !is_monomial(I) && dim(I) == 2
  @test associated_prime_heights(I) == [2, 4]
  pd = primary_decomposition(I^2)
  @test any(pair -> codim(pair[2]) == 3, pd)
  @test symbolic_power(I, 2; method=:primary) == I^2
  @test symbolic_power(I, 2) == I^2
  @test is_symbolic_equal_ordinary(I, 2)
  @test symbolic_power(I, 2; use_minimal_primes=true) == symbolic_power_join(J, 2)
  # This reproduces only the equality-based selection in Macaulay2's slow
  # fallback, not a Macaulay2 execution. It loses the new height-three prime.
  ass = [P for (_, P) in primary_decomposition(I)]
  equality_only = intersect([Q for (Q, P) in pd if any(A -> A == P, ass)]...)
  @test equality_only != I^2
end

# Hasse derivatives are computed directly from binomial coefficients, without
# elimination, primary decomposition, or SymbolicPowers internal helpers.
function reference_hasse_derivative(f, beta)
  R = parent(f)
  x = gens(R)
  result = zero(R)
  for (c, e) in zip(coefficients(f), exponents(f))
    all(e .>= beta) || continue
    factor = prod(binomial(e[i], beta[i]) for i in eachindex(e))
    result += c * factor * prod(x[i]^(e[i] - beta[i]) for i in eachindex(e))
  end
  return result
end

@testset "join against independent Hasse-derivative membership" begin
  for K in (QQ, GF(2), GF(3))
    R, (x, y, z) = polynomial_ring(K, [:x, :y, :z])
    # An invertible linear change of the radical triangle ideal.
    I = ideal(R, [(x + y) * y, (x + y) * z, y * z])
    for n in 2:3
      S = symbolic_power_join(I, n; check=false)
      orders = [
        collect(b) for b in Iterators.product(0:(n - 1), 0:(n - 1), 0:(n - 1)) if sum(b) < n
      ]
      differential_member(f) = all(b -> reference_hasse_derivative(f, b) in I, orders)
      @test all(differential_member, gens(S))
      for e in Iterators.product(0:3, 0:3, 0:3)
        sum(e) <= 4 || continue
        f = x^e[1] * y^e[2] * z^e[3]
        @test (f in S) == differential_member(f)
      end
    end
  end
  R, (x, y) = polynomial_ring(GF(2), [:x, :y])
  P = ideal(R, [x + y])
  @test !((x + y)^2 in symbolic_power_join(P, 3)) # ordinary derivatives would fail here
  @test symbolic_power_frobenius(P, 3) == P^3

  # F_2(s) is imperfect: x^2-s is prime, but the differential criterion fails.
  A, s = polynomial_ring(GF(2), :s)
  F = fraction_field(A)
  S, (u,) = polynomial_ring(F, [:u])
  Q = ideal(S, [u^2 - F(s)])
  @test_throws ArgumentError symbolic_power_join(Q, 2; check=false)
  @test symbolic_power_frobenius(Q, 2; check=false) == Q^2
end

@testset "squarefree products without constructing the ordinary power" begin
  R, (x, y, z, w, v) = polynomial_ring(QQ, [:x, :y, :z, :w, :v])
  pool = [x * y, z * w, x * z, y * w * v, x^2, y^2 * z]
  for mask in 1:63
    I = ideal(R, [pool[i] for i in eachindex(pool) if !iszero(mask & (1 << (i - 1)))])
    @test squarefree_in_codim(I) == squarefree_gens(I^codim(I))
  end
  I = ideal(R, [x * y, z * w, x * z, y * w * v])
  @test codim(I) == 2
  @test x * y * z * w in squarefree_in_codim(I)
  @test !(x * y * z * w * v in squarefree_in_codim(I)) # a redundant disjoint product
  @test squarefree_in_codim(ideal(R, elem_type(R)[])) == [one(R)]
  @test squarefree_in_codim(ideal(R, [one(R)])) == [one(R)]
  S, a = polynomial_ring(QQ, :a => 1:20)
  @test squarefree_in_codim(ideal(S, a)) == [prod(a)]
end

@testset "initial degree avoids unnecessary Groebner bases" begin
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  I = ideal(R, [x^3, y^2 + z^2, x * y * z])
  @test !Oscar.is_known(groebner_basis, I)
  @test initial_degree(I) == 2
  @test !Oscar.is_known(groebner_basis, I)
  U = ideal(R, [x, x + 1])
  @test initial_degree(U) == 0
  G, (u, v) = graded_polynomial_ring(QQ, [:u, :v]; weights=[2, 3])
  @test initial_degree(ideal(G, [u^3 + v^2, u^4])) == 6
end

@testset "known symbolic powers bypass backend exponent limits" begin
  R, _ = polynomial_ring(QQ, [:x, :y])
  Z, U = ideal(R, elem_type(R)[]), ideal(R, [one(R)])
  for I in (Z, U), method in SymbolicPowers._METHODS
    @test symbolic_power(I, typemax(Int); method) === I
  end
  @test symbolic_power_join(Z, typemax(Int)) === Z
  @test symbolic_power_join(U, typemax(Int)) === U
end

@testset "redundant decompositions of the power need no normalization" begin
  A, t = polynomial_ring(QQ, :t)
  K, a = number_field(t^2 + 1, "a")
  R, (x, y, z) = polynomial_ring(K, [:x, :y, :z])
  m = ideal(R, [x, y, z])
  examples = (
    ideal(R, [x * y, x * z, y * z]),
    intersect(ideal(R, [x - a * y, z]), ideal(R, [x, y])),
  )
  for I in examples, mp in (false, true)
    expected = symbolic_power(I, 2; method=:primary, use_minimal_primes=mp)
    power = I^2
    # Both added components contain I^2, since I ⊆ (x, y) ⊆ m. So the
    # decomposition stays valid; m is not an associated prime of I.
    redundant = vcat(
      primary_decomposition(power), [(m^2, m), (ideal(R, [x, y])^2, ideal(R, [x, y]))]
    )
    @test intersect(first.(redundant)...) == power
    set_attribute!(power, :primary_decomposition, redundant)
    @test SymbolicPowers._symbolic_power_primary(I, 2, mp, power) == expected
  end
end
