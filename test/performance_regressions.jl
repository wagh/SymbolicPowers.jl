@testset "reused structural data and exact monomial shortcuts" begin
  for K in (QQ, GF(5))
    R, (x, y, z) = polynomial_ring(K, [:x, :y, :z])
    I = ideal(R, [x * y, x * z, y * z])
    # Reuse the structure while changing exponents and symbolic conventions.
    for n in (2, 5, 3, 2), minimal in (false, true)
      S = symbolic_power(I, n; use_minimal_primes=minimal)
      @test S == intersect(ideal(R, [x, y])^n, ideal(R, [x, z])^n,
        ideal(R, [y, z])^n)
    end
    @test has_attribute(I, :symbolic_powers_minimal_primes)
    J = ideal(R, [x^3, x * y])
    for n in (1, 2, 4, 2), minimal in (false, true)
      expected = minimal ? ideal(R, [x^n]) : J^n
      @test symbolic_power(J, n; use_minimal_primes=minimal) == expected
      @test is_symbolic_equal_ordinary(J, n; use_minimal_primes=minimal) == !minimal
    end
    @test waldschmidt_constant(J) == 2
    @test waldschmidt_constant(J; use_minimal_primes=true) == 1
    # Independent quotient-module reference for the monomial defect shortcut.
    T, (a, b, c) = graded_polynomial_ring(K, [:a, :b, :c])
    examples = [ideal(T, [a * b, a * c, b * c]), ideal(T, [a^3, a * b]),
      ideal(T, [a^2, b^3]), ideal(T, [a^2 * b, a * b^2, c^2])]
    for H in examples, n in 1:4, minimal in (false, true)
      S = symbolic_power(H, n; use_minimal_primes=minimal)
      A, pr = quo(T, H^n)
      expected = length(minimal_generating_set(ideal(A, [pr(f) for f in gens(S)])))
      @test symbolic_defect(H, n; use_minimal_primes=minimal) == expected
    end
    # The squarefree H-description retains unused coordinates and agrees with
    # the convex hull of the known vertices plus the nonnegative orthant.
    V = matrix(QQ, [0 1 1; 1 0 1; 1 1 0; 1//2 1//2 1//2])
    @test symbolic_polyhedron(I) == convex_hull(V, identity_matrix(QQ, 3))
    @test waldschmidt_constant(I) == QQ(3, 2)
  end
  G, (a, b, c, unused) = graded_polynomial_ring(QQ, [:a, :b, :c, :unused];
    weights=[2, 3, 4, 7])
  @test waldschmidt_constant(ideal(G, [a * b, a * c, b * c])) == QQ(9, 2)
  @test waldschmidt_constant(ideal(G, [a^3, a * b])) == 5
  @test waldschmidt_constant(ideal(G, [a^3, a * b]); use_minimal_primes=true) == 2
  @test waldschmidt_constant(ideal(G, [a^3, b^2])) == 6
  @test symbolic_defect(ideal(G, [a * b, a * c, b * c]), 2) == 1
end
