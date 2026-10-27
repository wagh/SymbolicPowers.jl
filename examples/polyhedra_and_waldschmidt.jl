# Exact monomial polyhedra, certified upper bounds and resurgence lower bounds.
using Oscar, SymbolicPowers

R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
I = ideal(R, [x * y, x * z, y * z])
P = symbolic_polyhedron(I)
V = sort(Vector{QQFieldElem}.(collect(vertices(P))))
expected_vertices = sort([QQ.([1, 1, 0]), QQ.([1, 0, 1]), QQ.([0, 1, 1]),
  [QQ(1, 2), QQ(1, 2), QQ(1, 2)]])
@assert V == expected_vertices
@assert QQ.([1, 1, 1]) in 2 * P
@assert x * y * z in symbolic_power(I, 2)
println("Triangle symbolic polyhedron vertices: ", V)

exact = waldschmidt_constant(I)
upper1 = waldschmidt_upper_bound(I; sample_size=1)
upper2 = waldschmidt_upper_bound(I; sample_size=2)
@assert exact == QQ(3, 2)
@assert upper1 == 2
@assert upper2 == exact
println("Exact Waldschmidt constant: ", exact)
println("Sampled upper bounds (one, two powers): ", (upper1, upper2))

# N is primary and not squarefree. Its polyhedron contains a lattice point
# that is not an exponent of a monomial of N: polyhedra forget such holes.
N = ideal(R, [x^2, y^2])
@assert QQ.([1, 1, 0]) in symbolic_polyhedron(N)
@assert !(x * y in N)
@assert waldschmidt_constant(N) == 2
println("N = (x^2,y^2): (1,1,0) lies in its polyhedron, but xy is not in N.")

lower = lower_bound_resurgence(I; sample_size=5)
improved = lower_bound_resurgence(I; sample_size=5, use_waldschmidt=true)
@assert lower == QQ(6, 5)
@assert improved == QQ(4, 3)
println("Certified resurgence lower bounds: ", (lower, improved))

# For a nonmonomial homogeneous ideal, request a bound rather than an exact
# constant. Here linear complete-intersection theory separately proves 1.
Q = ideal(R, [x - y, y - z])
@assert initial_degree(Q) == 1
@assert waldschmidt_upper_bound(Q; sample_size=3) == 1
println("The sampled Waldschmidt upper bound for the point (1:1:1) is 1.")
