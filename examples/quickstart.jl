using Oscar, SymbolicPowers

R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
I = ideal(R, [x * y, x * z, y * z])
@assert symbolic_power(I, 2) == I^2 + ideal(R, [x * y * z])
@assert containment_problem(I, 2) == 3
@assert waldschmidt_constant(I) == QQ(3, 2)
@assert waldschmidt_upper_bound(I; sample_size=1) == 2
@assert symbolic_defect(I, 2) == 1

S, (a, b) = polynomial_ring(QQ, [:a, :b])
J = ideal(S, [a^2, a * b])
@assert symbolic_power(J, 2) == J^2
@assert symbolic_power(J, 2; use_minimal_primes=true) == ideal(S, [a^2])

F, (u, v) = polynomial_ring(GF(3), [:u, :v])
@assert symbolic_power_frobenius(ideal(F, [u, v]), 4) == ideal(F, [u, v])^4
