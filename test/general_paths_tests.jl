@testset "general nonmonomial paths" begin
  SP = SymbolicPowers

  # Reference: the general method on a fresh ideal object (no shared caches).
  reference(I, n, minimal) =
    symbolic_power(ideal(base_ring(I), gens(I)), n; method=:primary,
      use_minimal_primes=minimal)
  fresh(I) = ideal(base_ring(I), gens(I))
  # Unmixedness read off a primary decomposition, independent of the hull test.
  pd_unmixed(I) = begin
    pd = primary_decomposition(fresh(I))
    c = codim(I)
    all(((_, P),) -> codim(P) == c, pd)
  end

  @testset "unmixedness without primary decomposition" begin
    for K in (QQ, GF(7))
      R, (x, y, z, w) = polynomial_ring(K, [:x, :y, :z, :w])
      examples = [
        ideal(R, [x * z - y^2, x * w - y * z, y * w - z^2]),          # prime
        intersect(ideal(R, [x, y]), ideal(R, [z, w])),               # unmixed, radical
        intersect(ideal(R, [x + y + z]), ideal(R, [y - z, x + w])),  # mixed heights
        intersect(ideal(R, [x^2, x * y]), ideal(R, [z - 1])),        # embedded prime
        ideal(R, [x^2, y]),                                          # primary, nonreduced
        intersect(ideal(R, [x, y])^2, ideal(R, [z, w])),             # unmixed, nonreduced
        ideal(R, [x * (y - 1), x * z]),                              # (x) ∩ (y - 1, z)
      ]
      for I in examples
        @test SP._is_unmixed(fresh(I)) == pd_unmixed(I)
        J = fresh(I)
        @test big_height(J) ==
          maximum(codim(P) for (_, P) in primary_decomposition(fresh(I)))
        @test associated_prime_heights(fresh(I)) ==
          sort(unique([codim(P) for (_, P) in primary_decomposition(fresh(I))]))
      end
      # an already known decomposition is reused
      I = fresh(examples[3])
      primary_decomposition(I)
      @test !SP._is_unmixed(I)
      @test big_height(I) == 2
    end
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    @test big_height(ideal(R, [zero(R)])) == 0
    @test associated_prime_heights(ideal(R, [zero(R)])) == [0]
    @test_throws ArgumentError symbolic_power(intersect(ideal(R, [x]), ideal(R, [y, z])), 2;
      method=:unmixed)
  end

  @testset "automatic method choice" begin
    R, (x, y, z, w) = polynomial_ring(QQ, [:x, :y, :z, :w])
    twc = ideal(R, [x * z - y^2, x * w - y * z, y * w - z^2])
    mixed = intersect(ideal(R, [x^2 + y^2 - z]), ideal(R, [x - 1, y - z^2, w]))
    S, (a, b, c) = polynomial_ring(QQ, [:a, :b, :c])
    fermat = ideal(S, [a * (b^3 - c^3), b * (c^3 - a^3), c * (a^3 - b^3)])
    for minimal in (false, true)
      @test SP._choose_method(fresh(twc), minimal) == :unmixed
      @test SP._choose_method(fresh(mixed), minimal) == :components
    end
    # homogeneous of dimension one: the saturation for the default convention;
    # for the minimal-prime convention only through its shortcuts (the Fermat
    # points are not all rational), otherwise componentwise
    @test SP._choose_method(fresh(fermat), false) == :saturation
    @test SP._choose_method(fresh(fermat), true) == :components
    # minimal-prime convention through the saturation, with and without an
    # embedded irrelevant component
    pts = intersect(ideal(S, [a, b]), ideal(S, [b, c]), ideal(S, [a - c, b - c]))
    embedded = intersect(pts, ideal(S, [a, b, c])^3)
    for I in (fermat, pts, embedded)
      @test symbolic_power(fresh(I), 2; use_minimal_primes=true) == reference(I, 2, true)
    end
    @test symbolic_power(fresh(embedded), 2; use_minimal_primes=true) ==
      symbolic_power(fresh(pts), 2)
  end

  @testset "componentwise method agrees with the primary method" begin
    for K in (QQ, GF(5), GF(2, 2))
      R, (x, y, z, w) = polynomial_ring(K, [:x, :y, :z, :w])
      examples = [
        # radical, mixed heights (both targets isolated)
        intersect(ideal(R, [x + 2 * y + z]), ideal(R, [y - 2 * z, x + 3 * z])),
        intersect(ideal(R, [x^2 + y^2 - z]), ideal(R, [x - 1, y - z^2, w])),
        # three heights, one of them a point, and an unused direction
        intersect(
          ideal(R, [x * y - z^2]), ideal(R, [x, y]), ideal(R, [x - 1, y - 1, z - 1])
        ),
        # embedded prime (only nonisolated targets)
        ideal(R, [(x + z)^2, (x + z) * (y - z + 2 * x)]),
        # an embedded target together with an isolated one of the same height
        intersect(ideal(R, [x^2, x * y]), ideal(R, [z - 1])),
        # an embedded target, isolated targets of two heights
        intersect(ideal(R, [x^2, x * y]), ideal(R, [z - 1, w]), ideal(R, [y - 1])),
        # nonreduced isolated component next to a reduced one
        intersect(ideal(R, [x^2]), ideal(R, [y, z])),
        # (x) ∩ (y - 1, z), given by products
        ideal(R, [x * (y - 1), x * z]),
      ]
      for I in examples, minimal in (false, true), n in (K == QQ ? (1:3) : (2:2))
        S = SP._symbolic_power_components(fresh(I), n, minimal)
        @test S == reference(I, n, minimal)
        @test symbolic_power(fresh(I), n; use_minimal_primes=minimal) == S
      end
      for I in examples
        @test symbolic_power(fresh(I), 0) == ideal(R, [one(R)])
        @test symbolic_power(fresh(I), 1) == I
        @test symbolic_power(fresh(I), 1; use_minimal_primes=true) == minimal_part(fresh(I))
      end
      # one ideal object reused across exponents and conventions
      I = fresh(examples[6])
      for (n, minimal) in [(3, true), (2, false), (3, false), (2, true), (1, true)]
        @test symbolic_power(I, n; use_minimal_primes=minimal) == reference(I, n, minimal)
      end
    end
  end

  @testset "unmixed ideals over several fields and gradings" begin
    Qx, X = polynomial_ring(QQ, :X)
    Ki, i = number_field(X^2 + 1, :i)
    for K in (QQ, GF(7), Ki)
      R, (x, y, z, w) = polynomial_ring(K, [:x, :y, :z, :w])
      examples = [
        ideal(R, [x * z - y^2, x * w - y * z, y * w - z^2]),
        intersect(ideal(R, [x, y]), ideal(R, [z, w]), ideal(R, [x - z, y - w])),
        intersect(ideal(R, [x, y])^2, ideal(R, [z, w])),
        intersect(ideal(R, [x + 2 * y + z]), ideal(R, [y - 2 * z, x + 3 * z])),
        intersect(ideal(R, [x^2, x * y]), ideal(R, [z - 1])),
      ]
      for I in examples, minimal in (false, true)
        @test symbolic_power(fresh(I), 2; use_minimal_primes=minimal) ==
          reference(I, 2, minimal)
      end
    end
    G, (a, b, c, d) = graded_polynomial_ring(QQ, [:a, :b, :c, :d]; weights=[1, 2, 3, 4])
    for I in (ideal(G, [b^2 - a * c, a * d - b * c, c^2 - b * d]),
      intersect(ideal(G, [a^2 - b]), ideal(G, [b^2 - d, c - a^3])))
      for minimal in (false, true), n in 2:3
        @test symbolic_power(fresh(I), n; use_minimal_primes=minimal) ==
          reference(I, n, minimal)
      end
    end
  end

  @testset "mixed radical ideal at a larger exponent" begin
    # The primary method decomposes I^5 (about a minute with the default GTZ
    # algorithm). Independent value: (f) is prime and principal and Q is a prime
    # generated by a regular sequence, so their symbolic powers are ordinary
    # powers and I^(5) = (f^5) ∩ Q^5.
    R, (x, y, z, w) = polynomial_ring(QQ, [:x, :y, :z, :w])
    f = x^2 + y^2 - z
    Q = ideal(R, [x - 1, y - z^2, w])
    I = intersect(ideal(R, [f]), Q)
    @test symbolic_power(I, 5) == intersect(ideal(R, [f^5]), Q^5)
    @test symbolic_power(fresh(I), 5; use_minimal_primes=true) ==
      intersect(ideal(R, [f^5]), Q^5)
  end

  @testset "isolated components are hulled one at a time" begin
    # Curves with isolated (irrational or fat) points. Grouping the isolated
    # components by height made the points' part a long product generating set
    # (the hull of a zero-dimensional ideal is the ideal itself), and the final
    # intersection became very slow (line + 6 points, n = 3: over 15 minutes
    # for the primary method, 1.4 s per component).
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    pt(p) = ideal(R, [x - p[1], y - p[2], z - p[3]])
    examples = [
      intersect(ideal(R, [7y + 4z + 2, 7x - 5z + 8]),
        ideal(R, [z + 1, y^2 - 2, x - y]), pt([1, 1, 1])^2),
      intersect(ideal(R, [2y - z, 7x - 9z - 3]),
        ideal(R, [z - 4, y^2 - 3, x - 2y]), pt([1, -2, 5])^2),
      intersect(ideal(R, [x - y, z]), pt([1, 2, 3]), pt([2, -1, 1]),
        pt([0, 3, -1]), pt([-2, 1, 2])),
    ]
    for I in examples, mp in (false, true)
      @test symbolic_power(fresh(I), 2; use_minimal_primes=mp) ==
        symbolic_power(fresh(I), 2; method=:primary, use_minimal_primes=mp)
    end
  end
end
