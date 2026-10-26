using Test
using Oscar
using SymbolicPowers

# Most expected values below are taken from the test suite of the Macaulay2
# package SymbolicPowers (E. Grifo et al.) and from the examples of
# symbolicPowers.lib.

@testset "SymbolicPowers.jl" begin
  @testset "big height and associated primes" begin
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    @test big_height(ideal(R, [x, y^3, z^2])) == 3
    @test big_height(ideal(R, [x * (y^3 - z^3), y * (z^3 - x^3), z * (x^3 - y^3)])) == 2
    @test associated_prime_heights(ideal(R, [x^2, x * y])) == [1, 2]
    F, (a, b, c) = polynomial_ring(GF(2), [:a, :b, :c])
    @test big_height(ideal(F, [a, b])) == 2
    @test minimal_part(ideal(R, [x^2, x * y])) == ideal(R, [x])
    @test initial_degree(ideal(R, [x^3, y^2 + z^2, x * y * z])) == 2
  end

  @testset "symbolic powers (Macaulay2 tests)" begin
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    @test symbolic_power(ideal(R, [y - z, x + z]), 2) ==
      ideal(
      R, [y^2 - 2 * y * z + z^2, x * y - x * z + y * z - z^2, x^2 + 2 * x * z + z^2]
    )
    @test symbolic_power(ideal(R, [x]), 2) == ideal(R, [x^2])
    @test symbolic_power(ideal(R, [x + 1]), 2) == ideal(R, [x^2 + 2 * x + 1])
    @test symbolic_power(ideal(R, [x, y]), 2) == ideal(R, [x, y])^2
    @test symbolic_power(ideal(R, [x * y + x * z]), 2) == ideal(R, [(x * y + x * z)^2])

    S, (w, x, y, z) = polynomial_ring(QQ, [:w, :x, :y, :z])
    I = ideal(S, [x * y + 1, w * y * z])
    @test symbolic_power(I, 3) == ideal(
      S,
      [w^3 * z^3, w^2 * x * y * z^2 + w^2 * z^2,
        w * x^2 * y^2 * z + 2 * w * x * y * z + w * z,
        x^3 * y^3 + 3 * x^2 * y^2 + 3 * x * y + 1],
    )

    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    I = ideal(R, [x * (y^5 - z^5), y * (z^5 - x^5), z * (x^5 - y^5)])
    @test symbolic_power(I, 3; method=:complete_intersection) == symbolic_power(I, 3)
  end

  @testset "trivial cases" begin
    R, (x, y) = polynomial_ring(QQ, [:x, :y])
    I = ideal(R, [x * y])
    @test is_one(symbolic_power(I, 0))
    @test symbolic_power(I, 1) == I
    @test is_one(symbolic_power(ideal(R, [one(R)]), 3))
    @test_throws ArgumentError symbolic_power(I, -1)
    @test_throws ArgumentError symbolic_power(I, 2; method=:nonsense)
  end

  @testset "all methods agree" begin
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])

    # three coordinate points of P^2: squarefree monomial, dimension 1, radical
    I = ideal(R, [x * y, x * z, y * z])
    for n in 1:4
      ref = symbolic_power(I, n; method=:primary)
      @test symbolic_power(I, n) == ref
      @test symbolic_power(I, n; method=:monomial) == ref
      @test symbolic_power(I, n; method=:saturation) == ref
      @test symbolic_power(I, n; method=:unmixed) == ref
      @test symbolic_power(I, n; method=:join) == ref
    end
    @test symbolic_power(I, 2) == I^2 + ideal(R, [x * y * z])

    # the Fermat configuration of 12 points (examples of symbolicPowers.lib)
    J = ideal(R, [x * (y^3 - z^3), y * (z^3 - x^3), z * (x^3 - y^3)])
    ref = symbolic_power(J, 2; method=:primary)
    @test symbolic_power(J, 2) == ref
    @test symbolic_power(J, 2; method=:saturation) == ref
    @test symbolic_power(J, 2; method=:join) == ref

    # The monomial curve (t^3,t^4,t^5): a prime of height 2.
    P = ideal(R, [y^2 - x * z, x^2 * y - z^2, x^3 - y * z])
    ref = symbolic_power(P, 2; method=:primary)
    @test symbolic_power(P, 2) == ref
    @test symbolic_power(P, 2; method=:join) == ref
    @test symbolic_power(P, 2; method=:unmixed) == ref
    @test !(ref == P^2)

    # a non-squarefree monomial ideal with an embedded prime
    M = ideal(R, [x^3, x^2 * y^2, x * z^3, y^4, y^2 * z])
    for n in 2:3
      @test symbolic_power(M, n; method=:monomial) == symbolic_power(M, n; method=:primary)
      @test symbolic_power(M, n; method=:monomial, use_minimal_primes=true) ==
        symbolic_power(M, n; method=:primary, use_minimal_primes=true)
    end
  end

  @testset "graded rings" begin
    R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    I = ideal(R, [x * y, x * z, y * z])
    @test symbolic_power(I, 2) == I^2 + ideal(R, [x * y * z])
    @test symbolic_power(I, 3; method=:join) == symbolic_power(I, 3; method=:primary)
  end

  @testset "minimal primes versus associated primes" begin
    R, (x, y) = polynomial_ring(QQ, [:x, :y])
    I = ideal(R, [x^2, x * y])                 # (x) ∩ (x^2, y)
    @test symbolic_power(I, 2) == I^2        # (x, y) is associated
    @test symbolic_power(I, 2; use_minimal_primes=true) == ideal(R, [x^2])
  end

  @testset "join" begin
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    @test symbolic_join(ideal(R, [x, y]), ideal(R, [x, z])) == ideal(R, [x])
    @test symbolic_power_join(ideal(R, [x, y]), 3) == ideal(R, [x, y])^3
    @test_throws ArgumentError symbolic_power(ideal(R, [x^2]), 2; method=:join)
  end

  @testset "positive characteristic" begin
    # Singular's primary decomposition can fail over very small prime fields,
    # so the Frobenius method is compared with the join method there.
    for p in (2, 3, 5)
      F, (x, y, z) = polynomial_ring(GF(p), [:x, :y, :z])
      P = ideal(F, [y^2 - x * z, x^2 * y - z^2, x^3 - y * z])
      for n in 2:4
        @test symbolic_power(P, n; method=:frobenius) == symbolic_power(P, n; method=:join)
      end
      L = ideal(F, [x, y])
      for n in 1:5
        @test symbolic_power_frobenius(L, n) == L^n
      end
    end
    # over a larger prime field the general method is available as a reference
    F, (x, y, z) = polynomial_ring(GF(32003), [:x, :y, :z])
    P = ideal(F, [y^2 - x * z, x^2 * y - z^2, x^3 - y * z])
    for n in 2:3
      @test symbolic_power(P, n; method=:join) == symbolic_power(P, n; method=:primary)
    end
    R, (x, y) = polynomial_ring(QQ, [:x, :y])
    @test_throws ArgumentError symbolic_power_frobenius(ideal(R, [x]), 2)
  end

  @testset "containment problem" begin
    S, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    @test is_symbolic_power_contained_in_power(ideal(S, [x, y, z]), 2, 2)

    R, (x, y) = polynomial_ring(QQ, [:x, :y])
    I = ideal(R, [x])
    @test !is_symbolic_power_contained_in_power(I, 2, 3)
    @test is_symbolic_power_contained_in_power(I, 2, 2)
    @test is_symbolic_power_contained_in_power(I, 3, 2)

    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    @test containment_problem(ideal(R, [x * y, x * z, y * z]), 2) == 3
    J = ideal(R, [x * (y^3 - z^3), y * (z^3 - x^3), z * (x^3 - y^3)])
    @test containment_problem(J, 2) == 4
    @test containment_problem(J, 4; in_symbolic=true) == 2
    @test lower_bound_resurgence(ideal(R, [x * y, x * z, y * z])) == QQ(6, 5)
  end

  @testset "symbolic equals ordinary" begin
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    @test !is_symbolic_equal_ordinary(ideal(R, [x * y, x * z, y * z]), 2)
    @test is_symbolic_equal_ordinary(ideal(R, [x * z, y * z]), 2)
    F, (a, b) = polynomial_ring(GF(3), [:a, :b])
    @test is_symbolic_equal_ordinary(ideal(F, [a]), 3)
  end

  @testset "packing problem" begin
    R5, (w, x, y, z) = polynomial_ring(GF(5), [:w, :x, :y, :z])
    @test squarefree_gens(ideal(R5, [y^2 * z, x * y * w, z * w^3])) == [w * x * y]

    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    @test isempty(squarefree_gens(ideal(R, [x^2 * z, x * y^8, z^3])))
    @test isempty(squarefree_in_codim(ideal(R, [x, y^2])))
    F, (a, b, c) = polynomial_ring(GF(2), [:a, :b, :c])
    @test squarefree_in_codim(ideal(F, [a, b])) == [a * b]

    R17, _ = polynomial_ring(GF(17), [:x])
    @test is_konig(ideal(R17, [one(R17)]))
    T, _ = polynomial_ring(QQ, [:y, :z])
    @test is_konig(ideal(T, elem_type(T)[]))

    I = ideal(R, [x * y, x * z, y * z])
    @test !is_konig(I)
    @test !is_packed(I)
    nps = non_packed_substitutions(I)
    @test length(nps) == 1 && isempty(nps[1].zero) && isempty(nps[1].one)

    S, (x, y, z, a, b) = polynomial_ring(QQ, [:x, :y, :z, :a, :b])
    J = intersect(ideal(S, [x, y]), ideal(S, [x, z]), ideal(S, [z, a]),
      ideal(S, [y, a]), ideal(S, [x, b]))
    @test is_packed(J)
    @test isempty(non_packed_substitutions(J))

    U, v = polynomial_ring(QQ, :x => 1:6)
    K = intersect(ideal(U, [v[1], v[2]]), ideal(U, [v[2], v[3]]), ideal(U, [v[3], v[1]]),
      ideal(U, [v[3], v[4]]), ideal(U, [v[4], v[5]]), ideal(U, [v[5], v[6]]),
      ideal(U, [v[6], v[4]]))
    witnesses = non_packed_substitutions(K)
    @test !is_packed(K)
    # the trivial substitution is reported exactly when K itself is not König
    @test any(s -> isempty(s.zero) && isempty(s.one), witnesses) == !is_konig(K)
    # the four substitutions found by Macaulay2's noPackedAllSubs
    m2 = [([1], [2, 3]), ([2], [1, 3]), ([5], [4, 6]), ([6], [4, 5])]
    found = [
      (sort([findfirst(==(t), v) for t in s.zero]),
        sort([findfirst(==(t), v) for t in s.one])) for s in witnesses
    ]
    @test all(in(found), m2)
    @test length(non_packed_substitutions(K; first_only=true)) == 1
  end

  @testset "symbolic polyhedron and Waldschmidt constant" begin
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    I = ideal(R, [x * y, y * z, x * z])
    V = sort(Vector{QQFieldElem}.(collect(vertices(symbolic_polyhedron(I)))))
    @test V == sort([QQ.([1, 1, 0]), QQ.([1, 0, 1]), QQ.([0, 1, 1]),
      [QQ(1, 2), QQ(1, 2), QQ(1, 2)]])
    @test waldschmidt_constant(I) == QQ(3, 2)
    # the sampled bound for a non-monomial ideal, compared with the exact value
    J = ideal(R, [x - y, y - z]) # a point; Waldschmidt constant 1
    @test waldschmidt_upper_bound(J) == 1
    @test_throws ArgumentError waldschmidt_constant(J)

    # monomial curve (t^9, t^11, t^13) in a weighted polynomial ring
    G, (a, b, c) = graded_polynomial_ring(QQ, [:a, :b, :c]; weights=[9, 11, 13])
    T, (t,) = polynomial_ring(QQ, [:t])
    phi = hom(G, T, [t^9, t^11, t^13])
    K = kernel(phi)
    @test waldschmidt_upper_bound(K) == 22 # Finite upper bound, not an exact-value claim.
    L = leading_ideal(K; ordering=wdegrevlex(G, [9, 11, 13]))
    @test waldschmidt_constant(L) == 22
  end

  @testset "symbolic defect and asymptotic regularity" begin
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    I = ideal(R, [x * y, x * z, y * z])
    @test symbolic_defect(I, 1) == 0
    @test symbolic_defect(I, 2) == 1
    @test min_degree_symbolic_power(I, 2) == 3
    @test_throws ArgumentError asymptotic_regularity(ideal(R, [x, y]); sample_size=3)
    @test asymptotic_regularity(ideal(R, [x, y]); sample_size=3, approximate=true) == 1
    @test symbolic_regularity_samples(ideal(R, [x, y]); sample_size=3) ==
      [(power=n, regularity=n, ratio=QQ(1)) for n in 1:3]
  end

  include("regressions.jl")
  include("monomial_oracle.jl")
  include("review_regressions.jl")
  include("performance_regressions.jl")
  include("monomial_engine_tests.jl")
end
