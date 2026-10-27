@testset "large exponents preserve exact arithmetic" begin
  SP = SymbolicPowers
  h = div(typemax(Int), 2) + 1
  # Total degrees can exceed Int even though every exponent fits in Int.
  @test SP._mono_degree([h, h]) == 2 * ZZ(h)
  @test SP._mono_minimalize([[h, h], [1, 0]]) == [[1, 0]]
  A = SP._mono_minimalize([[h, 0], [1, 1]])
  @test SP._mono_intersect(A, [[0, h]]) == [[1, h]]
  @test SP._mono_intersect_prime_power([[h, h]], [1, 2], 1) == [[h, h]]
  # Mask thresholds must remain positive and increasing near typemax(Int).
  mk = SP._MonoMasker(1, [[typemax(Int)]])
  @test issorted(mk.thresholds[1]) && all(t -> 1 <= t <= typemax(Int), mk.thresholds[1])
  @test SP._mono_in([typemax(Int)], [[h]])
  @test_throws OverflowError SP._mono_product([[typemax(Int)]], [[1]])

  R, (x, y) = polynomial_ring(QQ, [:x, :y])
  I = ideal(R, [x, x^h * y^h]) # The second generator is redundant: I = (x).
  @test SP._minimal_monomial_generators(I) == [x]
  @test is_packed(I)
  @test symbolic_power(I, 3) == ideal(R, [x^3])

  # A wrapped fat-point exponent formerly returned P^4 for n*4 = 2^wordsize+4.
  n = div(typemax(Int), 2) + 2
  @test 4 * ZZ(n) > typemax(Int)
  for K in (QQ, GF(7))
    T, (a, b, c) = polynomial_ring(K, [:a, :b, :c])
    P = ideal(T, [a - b, b - c])
    J = P^4
    comps = SP._rational_point_components(J)
    @test comps !== nothing && all(pair -> pair[2] == 4, comps)
    @test symbolic_power(J, 2) == P^8
    for mp in (false, true), method in (:auto, :saturation)
      @test_throws OverflowError symbolic_power(J, n; method, use_minimal_primes=mp)
    end
  end
end

@testset "monomial unit checks do not require a Groebner basis" begin
  for K in (QQ, GF(7))
    R, (x, y) = polynomial_ring(K, [:x, :y])
    for (I, expected) in ((ideal(R, [x^2, y^3]), false),
      (ideal(R, [zero(R)]), false), (ideal(R, [x, one(R)]), true))
      @test !Oscar.is_known(groebner_basis, I)
      @test SymbolicPowers._is_unit_ideal(I) == expected
      @test !Oscar.is_known(groebner_basis, I)
    end
    @test SymbolicPowers._is_unit_ideal(ideal(R, [x, x + 1]))
  end
end

@testset "Frobenius complete-intersection certificate" begin
  for p in (3, 101)
    R, (x, y, z, w) = polynomial_ring(GF(p), [:x, :y, :z, :w])
    # Quotient rings are polynomial rings or a quadric domain, hence prime.
    # These include nonlinear complete intersections as well as linear ones.
    primes = [ideal(R, [x + y, z]), ideal(R, [x + y, z^2 + w]),
      ideal(R, [x * y - z^2, w])]
    for P in primes, n in 2:3
      @test symbolic_power_frobenius(P, n) == P^n
    end
    @test_throws ArgumentError symbolic_power_frobenius(primes[1]^2, 2)
  end
end

@testset "polynomial rings with no variables" begin
  for K in (QQ, GF(3))
    R, _ = polynomial_ring(K, Symbol[])
    Z, U = ideal(R, elem_type(R)[]), ideal(R, [one(R)])
    for I in (Z, U), J in (Z, U)
      S = symbolic_join(I, J)
      @test base_ring(S) === R
      @test is_zero(S) == (is_zero(I) && is_zero(J))
      @test is_one(S) == (is_one(I) || is_one(J))
    end
    for I in (Z, U), n in 0:2
      S = symbolic_power_join(I, n)
      @test n == 0 ? is_one(S) : S === I
      @test n == 0 ? is_one(symbolic_power(I, n)) : symbolic_power(I, n) === I
    end
    @test big_height(Z) == 0
    @test initial_degree(U) == 0
    @test is_packed(Z) && is_packed(U)
    if characteristic(K) > 0
      for n in 0:2
        S = symbolic_power_frobenius(Z, n)
        @test n == 0 ? is_one(S) : S === Z
      end
      @test_throws ArgumentError symbolic_power_frobenius(U, 2)
    end
  end
end

@testset "principal monomial powers at the Int boundary" begin
  R, (x,) = polynomial_ring(QQ, [:x])
  I = ideal(R, [x])
  n = typemax(Int)
  @test SymbolicPowers._mono_power([zeros(Int, 2)], n, 2) == [zeros(Int, 2)]
  @test is_symbolic_equal_ordinary(I, n)
  @test containment_problem(I, n; in_symbolic=true, use_minimal_primes=true) == n
  @test containment_problem(I, n; use_minimal_primes=true, max_power=n) == n
end
