# Number fields, finite extension fields, weighted degrees and multigradings.
using Oscar, SymbolicPowers

A, t = polynomial_ring(QQ, :t)
K, s = number_field(t^2 - 2, :s)
R, (x, y, z) = polynomial_ring(K, [:x, :y, :z])
P = ideal(R, [x - s * y, z])
@assert symbolic_power(P, 2; method=:primary) == P^2
@assert symbolic_power_join(P, 2) == P^2
println("Over QQ(s), s^2 = 2: the prime (x-s*y,z) has P^(2) = P^2.")

F, a = finite_field(3, 2, "a")
T, (u, v, w) = polynomial_ring(F, [:u, :v, :w])
Q = ideal(T, [u - a * v, w])
@assert symbolic_power_frobenius(Q, 2) == Q^2
@assert symbolic_power_join(Q, 2) == Q^2
@assert frobenius_power(Q, 3) == ideal(T, [u^3 - a^3 * v^3, w^3])
println("Over GF(9), Frobenius also raises coefficients: (u-a*v)^3 = u^3-a^3*v^3.")

G, (b, c) = graded_polynomial_ring(QQ, [:b, :c]; weights=[2, 3])
I = ideal(G, [b^3, c^2])
@assert initial_degree(I) == 6
@assert min_degree_symbolic_power(I, 2) == 12
@assert waldschmidt_constant(I) == 6
println("Weights (2,3): alpha((b^3,c^2)) = 6 and its Waldschmidt constant is 6.")

# Symbolic powers also make sense in a multigrading. The degree-valued APIs
# require a positive Z-grading, so no scalar initial degree is requested here.
H, (h1, h2, h3) = graded_polynomial_ring(QQ, [:h1, :h2, :h3];
  weights=[[1, 0], [1, 0], [0, 1]])
L = ideal(H, [h1 + h2])
@assert symbolic_power(L, 2) == L^2
println("A Z^2-graded principal prime still has L^(2) = L^2.")
