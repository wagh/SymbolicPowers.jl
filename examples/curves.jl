# An explicit symbolic-square element for the monomial curve (t^3,t^4,t^5).
# Weighted grading makes its binomial equations homogeneous.
using Oscar, SymbolicPowers

R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]; weights=[3, 4, 5])
A = y^2 - x * z
B = x^2 * y - z^2
C = x^3 - y * z
P = ideal(R, [A, B, C])
f = x^5 + x * y^3 + z^3 - 3 * x^2 * y * z

# x*f is in P^2, while x is outside the prime P. Thus f is in P^(2).
# Its weighted degree is 15, below the least degree 16 of P^2.
@assert x * f == C^2 + A * B
@assert !(x in P)
@assert initial_degree(P) == 8
S = symbolic_power(P, 2; method=:unmixed)
@assert S == symbolic_power(P, 2; method=:primary)
@assert f in S
@assert !(f in P^2)
@assert min_degree_symbolic_power(P, 2) == 15
println("Monomial curve: alpha(P) = 8; alpha(P^(2)) = 15 < alpha(P^2) = 16.")
println("An explicit element of P^(2) outside P^2: ", f)

# The twisted cubic gives a second family, in four variables and the standard
# grading. Compare methods by ideal equality, not by generator-list equality.
T, (a, b, c, d) = graded_polynomial_ring(QQ, [:a, :b, :c, :d])
J = ideal(T, [a * c - b^2, a * d - b * c, b * d - c^2])
reference = symbolic_power(J, 2; method=:primary)
@assert symbolic_power(J, 2) == reference
@assert symbolic_power(J, 2; method=:unmixed) == reference
println("Twisted cubic: automatic, primary and unmixed methods agree.")
