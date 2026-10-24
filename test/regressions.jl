@testset "zero ideals and exponent boundaries" begin
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  Z, U = ideal(R, elem_type(R)[]), ideal(R, [one(R)])
  for m in 0:3, n in 0:3
    @test is_symbolic_power_contained_in_power(Z, m, n) == (n == 0 || m > 0)
    @test is_symbolic_power_contained_in_power(U, m, n)
  end
  @test is_one(symbolic_power(Z, 0))
  @test is_zero(symbolic_power(Z, 1))
  @test is_zero(symbolic_power(Z, 3))
  @test associated_prime_heights(Z) == [0]
  @test big_height(Z) == 0
  @test minimal_part(Z) == Z
  @test_throws ArgumentError big_height(U)
  @test_throws ArgumentError initial_degree(Z)
  @test_throws ArgumentError waldschmidt_constant(U)
  @test_throws ArgumentError waldschmidt_constant(Z)
  for I in (Z, U)
    @test symbolic_defect(I, 0) == 0
    @test symbolic_defect(I, 2) == 0
    @test_throws ArgumentError containment_problem(I, 1)
    @test_throws ArgumentError lower_bound_resurgence(I)
    @test_throws ArgumentError symbolic_regularity_samples(I)
  end
  @test_throws ArgumentError containment_problem(ideal(R, [x]), 2; max_power=1)
  @test_throws ArgumentError containment_problem(ideal(R, [x]), 1; max_power=0)
  @test_throws ArgumentError lower_bound_resurgence(ideal(R, [x]); sample_size=0)
  @test_throws ArgumentError waldschmidt_upper_bound(ideal(R, [x]); sample_size=0)
end

@testset "embedded components and non-equidimensional support" begin
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  I = ideal(R, [x^2, x * y])
  @test symbolic_power(I, 1; use_minimal_primes=true) == ideal(R, [x])
  for n in 1:3
    @test symbolic_power(I, n; method=:primary) == I^n
  end
  # In three variables I has dimension 2; saturation must refuse it.
  @test_throws ArgumentError symbolic_power(I, 2; method=:saturation)

  S, (a, b) = polynomial_ring(QQ, [:a, :b])
  J = ideal(S, [a^2, a * b])
  for n in 1:3
    @test symbolic_power(J, n; method=:saturation) == J^n
    @test symbolic_power(J, n; method=:saturation, use_minimal_primes=true) ==
      ideal(S, [a^n])
  end
  @test containment_problem(J, 2; use_minimal_primes=true, max_power=4) == 4
  @test_throws ArgumentError containment_problem(J, 2; use_minimal_primes=true, max_power=3)
  @test containment_problem(J, 1; in_symbolic=true, use_minimal_primes=true) == 0
  @test symbolic_defect(J, 1; use_minimal_primes=true) == 1

  # Incomparable minimal primes of different heights must both survive.
  H = ideal(R, [x * y, x * z])
  for n in 1:3
    expected = intersect(ideal(R, [x^n]), ideal(R, [y, z])^n)
    @test symbolic_power(H, n; method=:primary) == expected
    @test symbolic_power(H, n; method=:monomial) == expected
    @test symbolic_power(H, n; method=:join) == expected
  end
  @test_throws ArgumentError symbolic_power(H, 2; method=:unmixed)
  @test_throws ArgumentError symbolic_power(H, 2; method=:complete_intersection)
end

@testset "checked join and Frobenius contracts" begin
  R, (x, y) = polynomial_ring(QQ, [:x, :y])
  I = ideal(R, [x^2])
  @test symbolic_join(I, ideal(R, [x, y])^2) == ideal(R, [x^3])
  @test symbolic_power(I, 2; method=:primary) == ideal(R, [x^4])
  @test_throws ArgumentError symbolic_power_join(I, 2)
  @test symbolic_power_join(I, 2; check=false) == ideal(R, [x^3])
  S, (a, b) = polynomial_ring(QQ, [:a, :b]; cached=false)
  @test_throws ArgumentError symbolic_join(ideal(R, [x]), ideal(S, [a]))
  @test_throws ArgumentError frobenius_power(I, 2)
  @test frobenius_power(I, 1) === I
  for p in (2, 3, 5)
    F, (a, b) = polynomial_ring(GF(p), [:a, :b])
    J = ideal(F, [a, b])
    redundant = ideal(F, [a, b, a + b])
    @test frobenius_power(J, p) == frobenius_power(redundant, p)
    @test frobenius_power(J, p^2) == ideal(F, [a^(p^2), b^(p^2)])
    @test_throws ArgumentError frobenius_power(J, p + 1)
    @test_throws ArgumentError symbolic_power_frobenius(J^2, 2)
    for n in (p - 1, p, p + 1)
      @test symbolic_power_frobenius(ideal(F, [a]), n) == ideal(F, [a^n])
      @test symbolic_power_frobenius(J, n) == J^n
    end
  end
end

@testset "König supports, packing domain and lazy enumeration" begin
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  @test is_konig(ideal(R, [x^2, y^2]))
  @test is_konig(ideal(R, [x^2 * y, z^3]))
  @test !is_konig(ideal(R, [x^2 * y, x * z^2, y * z]))
  @test_throws ArgumentError is_konig(ideal(R, [x + y]))
  @test_throws ArgumentError is_packed(ideal(R, [x^2, y^2]))
  @test_throws ArgumentError non_packed_substitutions(ideal(R, [x^2]))
  # An eager list of 3^40 pairs used to overflow/allocate before finding a witness.
  W, v = polynomial_ring(QQ, :v => 1:40)
  T = ideal(W, [v[1] * v[2], v[1] * v[3], v[2] * v[3]])
  @test !is_packed(T)
  @test length(non_packed_substitutions(T; first_only=true)) == 1
  # only the three active variables are substituted, so the full list is small
  all_w = non_packed_substitutions(T)
  @test length(all_w) == 1 && isempty(all_w[1].zero) && isempty(all_w[1].one)
  @test is_packed(ideal(W, [v[1] * v[2], v[3] * v[4]]))
  @test is_packed(ideal(W, [one(W)])) && is_packed(ideal(W, elem_type(W)[]))
end

@testset "grading and exact versus sampled invariants" begin
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  I = ideal(R, [x^2, y^2])
  @test !(x * y in I)
  @test QQ.([1, 1, 0]) in symbolic_polyhedron(I)
  @test waldschmidt_constant(I) == 2
  T = ideal(R, [x * y, x * z, y * z])
  @test waldschmidt_constant(T) == QQ(3, 2)
  @test waldschmidt_upper_bound(T; sample_size=1) == 2
  @test waldschmidt_upper_bound(T; sample_size=2) == QQ(3, 2)
  @test symbolic_defect(T, 0) == 0
  @test_throws ArgumentError initial_degree(ideal(R, [x + 1]))
  @test_throws ArgumentError waldschmidt_upper_bound(ideal(R, [x + 1]))
  @test initial_degree(ideal(R, [x, x + y^2])) == 1

  G, (a, b) = graded_polynomial_ring(QQ, [:a, :b]; weights=[2, 3])
  @test initial_degree(ideal(G, [a^3, b^2])) == 6
  @test waldschmidt_constant(ideal(G, [a^3, b^2])) == 6
  W, (u, v) = graded_polynomial_ring(QQ, [:u, :v]; weights=[typemax(Int), 1])
  @test_throws Union{OverflowError,InexactError} initial_degree(ideal(W, [u^2]))
  @test_throws OverflowError SymbolicPowers._weighted_degree(x^2, [typemax(Int), 1, 1])
  @test_throws ArgumentError symbolic_regularity_samples(ideal(G, [a, b]))
  N, (a, b) = graded_polynomial_ring(QQ, [:a, :b]; weights=[0, 1])
  @test_throws ArgumentError initial_degree(ideal(N, [a]))
  M, (a, b, c) = graded_polynomial_ring(QQ, [:a, :b, :c];
    weights=[[1, 0], [1, 0], [0, 1]])
  @test_throws ArgumentError initial_degree(ideal(M, [a, b]))
  @test symbolic_power(ideal(M, [a + b]), 2) == ideal(M, [(a + b)^2])
  A, _ = polynomial_ring(ZZ, [:a, :b])
  @test_throws ArgumentError symbolic_power(ideal(A, gens(A)), 2)

  # Corona of K_3 with two pendant edges at each vertex: the symbolic
  # regularity slope is 9/2 although the first ratio is 4 (DNS, Ex. 3.16).
  C, (a, b, c, p, q, r, s, t, u) = polynomial_ring(QQ,
    [:a, :b, :c, :p, :q, :r, :s, :t, :u])
  J = ideal(C, [a * b * c, b * c * p * q, a * c * r * s, a * b * t * u])
  data = symbolic_regularity_samples(J; sample_size=3)
  @test [d.regularity for d in data] == [4, 9, 13]
  @test minimum(d.ratio for d in data) == 4 < QQ(9, 2)
end

@testset "nonrational residue points and finite extension fields" begin
  R, (x, y) = polynomial_ring(QQ, [:x, :y])
  P = ideal(R, [x^2 + 1, y]) # A maximal ideal with residue field QQ(i).
  for method in (:primary, :unmixed, :complete_intersection, :join)
    @test symbolic_power(P, 2; method) == P^2
  end
  K, a = finite_field(2, 2, "a")
  F, (u, v) = polynomial_ring(K, [:u, :v])
  Q = ideal(F, [u + a * v])
  @test symbolic_power_join(Q, 3) == Q^3
  @test symbolic_power_frobenius(Q, 3) == Q^3
  @test frobenius_power(ideal(F, [u, v]), 4) == ideal(F, [u^4, v^4])
  @test_throws ArgumentError symbolic_power_frobenius(Q, typemax(Int))
end

# Brute-force references built only from OSCAR ideal operations.
function reference_is_konig(I)
  (is_one(I) || is_zero(I)) && return true
  c = codim(I)
  G = gens(groebner_basis(I; complete_reduction=true))   # minimal monomial generators
  coprime(a, b) = is_one(gcd(a, b))
  for S in Iterators.product(ntuple(_ -> eachindex(G), c)...)
    (issorted(S) && length(unique(S)) == c) || continue # strictly increasing indices
    all(coprime(G[S[i]], G[S[j]]) for i in 1:c for j in (i + 1):c) && return true
  end
  return false
end

function reference_minor(I, Z, O)
  R = base_ring(I)
  imgs = [i in Z ? zero(R) : (i in O ? one(R) : x) for (i, x) in enumerate(gens(R))]
  g = filter(!iszero, [evaluate(f, imgs) for f in gens(I)])
  return isempty(g) ? ideal(R, elem_type(R)[]) : ideal(R, g)
end

@testset "combinatorial König and packing against ideal computations" begin
  R, x = polynomial_ring(QQ, :x => 1:4)
  pool = [
    x[1] * x[2], x[2] * x[3], x[3] * x[4], x[1] * x[4], x[1] * x[3], x[2] * x[3] * x[4]
  ]
  nonsquarefree = [x[1]^2, x[2]^3 * x[3], x[3] * x[4]^2, x[1] * x[4]]
  for mask in 1:(2^length(pool) - 1)
    I = ideal(R, [pool[i] for i in eachindex(pool) if !iszero(mask & (1 << (i - 1)))])
    @test is_konig(I) == reference_is_konig(I)
    # packing over all substitutions of the active variables, i.e. those
    # occurring in a minimal generator (a redundant generator such as x2*x3*x4
    # next to x2*x3 does not make x4 active)
    mingens = gens(groebner_basis(I; complete_reduction=true))
    active = findall(i -> any(f -> degree(f, i) > 0, mingens), 1:4)
    expected = Tuple{Vector{Int},Vector{Int}}[]
    for code in Iterators.product(ntuple(_ -> 0:2, length(active))...)
      Z = active[findall(==(1), collect(code))]
      O = active[findall(==(2), collect(code))]
      reference_is_konig(reference_minor(I, Z, O)) || push!(expected, (Z, O))
    end
    found = [
      (
        sort([findfirst(==(t), x) for t in w.zero]),
        sort([findfirst(==(t), x) for t in w.one]),
      )
      for w in non_packed_substitutions(I)
    ]
    @test sort(found) == sort(expected)
    @test is_packed(I) == isempty(expected)
  end
  for mask in 1:(2^length(nonsquarefree) - 1)
    I = ideal(
      R,
      [
        nonsquarefree[i] for i in eachindex(nonsquarefree) if !iszero(mask & (1 << (i - 1)))
      ],
    )
    @test is_konig(I) == reference_is_konig(I)
  end
end

@testset "containment searches against a linear scan" begin
  R, (x, y, z, w) = polynomial_ring(QQ, [:x, :y, :z, :w])
  linear(I, n; kw...) = first(m for m in n:typemax(Int)
        if is_symbolic_power_contained_in_power(I, m, n; kw...))
  largest(I, n; kw...) = last(k for k in 0:n
                              if is_subset(symbolic_power(I, n; kw...), I^k))
  examples = [
    (ideal(R, [x * y, x * z, y * z]), 1:4),                                 # big height 2
    (intersect(ideal(R, [x * y, x * z, y * z]), ideal(R, [x, y, z, w])^3), 1:3), # embedded maximal ideal, big height 4
    (ideal(R, [x * (y^3 - z^3), y * (z^3 - x^3), z * (x^3 - y^3)]), 1:2),   # Fermat points
    (ideal(R, [x^2, x * y]), 1:3),
  ]
  for (I, ns) in examples, n in ns
    @test containment_problem(I, n) == linear(I, n)
    @test containment_problem(I, n; in_symbolic=true) == largest(I, n)
  end
  J = ideal(R, [x^2, x * y])
  for n in 1:2
    @test containment_problem(J, n; use_minimal_primes=true, max_power=8) ==
      linear(J, n; use_minimal_primes=true)
    @test containment_problem(J, n; in_symbolic=true, use_minimal_primes=true) ==
      largest(J, n; use_minimal_primes=true)
  end
  # the resurgence bound shares symbolic powers between the searches
  T = ideal(R, [x * y, x * z, y * z])
  @test lower_bound_resurgence(T; sample_size=5) == QQ(6, 5)
  @test lower_bound_resurgence(T; sample_size=3, use_waldschmidt=true) == QQ(4, 3)
end

@testset "symbolic equals ordinary, both conventions" begin
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  I = ideal(R, [x^2, x * y])                      # (x) ∩ (x^2, y)
  @test is_symbolic_equal_ordinary(I, 1)
  @test !is_symbolic_equal_ordinary(I, 1; use_minimal_primes=true)
  @test is_symbolic_equal_ordinary(I, 2)
  for n in 0:3,
    J in (I, ideal(R, [x * y, x * z, y * z]), ideal(R, [x * z, y * z]),
      ideal(R, [y^2 - x * z, x^2 * y - z^2, x^3 - y * z]))

    for mp in (false, true)
      @test is_symbolic_equal_ordinary(J, n; use_minimal_primes=mp) ==
        (symbolic_power(J, n; method=:primary, use_minimal_primes=mp) == J^n)
    end
  end
end

@testset "homogeneity is a property of the ideal, not of its generators" begin
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  I = ideal(R, [x * y, x * z, y * z + x * y * z])  # = (xy, xz, yz)
  T = ideal(R, [x * y, x * z, y * z])
  @test symbolic_power(I, 2; method=:saturation) == symbolic_power(T, 2)
  @test SymbolicPowers._choose_method(
    ideal(R, [x * (y^3 - z^3), y * (z^3 - x^3),
      z * (x^3 - y^3) + x * (y^3 - z^3) * z]), false) == :saturation
  @test_throws ArgumentError symbolic_power(ideal(R, [x - 1, y]), 2; method=:saturation)
end
