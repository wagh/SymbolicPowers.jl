# Associated-prime versus minimal-prime symbolic powers, and mixed heights.
using Oscar, SymbolicPowers

R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
I = ideal(R, [x^2, x * y]) # (x) intersect (x^2,y); (x,y) is embedded.
@assert associated_prime_heights(I) == [1, 2]
@assert big_height(I) == 2
@assert minimal_part(I) == ideal(R, [x])

for n in 1:3
  @assert symbolic_power(I, n) == I^n
  @assert symbolic_power(I, n; use_minimal_primes=true) == ideal(R, [x^n])
end
@assert is_symbolic_equal_ordinary(I, 1)
@assert !is_symbolic_equal_ordinary(I, 1; use_minimal_primes=true)
println("Embedded example: heights = ", associated_prime_heights(I))
println("Default symbolic square: ", symbolic_power(I, 2))
println("Minimal-prime symbolic square: ", symbolic_power(I, 2; use_minimal_primes=true))

# A radical ideal may have minimal primes of different heights. Both must
# survive: taking only the largest-dimensional part would lose a component.
J = ideal(R, [x * y, x * z]) # (x) intersect (y,z).
expected = intersect(ideal(R, [x^2]), ideal(R, [y, z])^2)
@assert associated_prime_heights(J) == [1, 2]
@assert minimal_part(J) == J
@assert symbolic_power(J, 2) == expected
@assert symbolic_power(J, 2; use_minimal_primes=true) == expected
println("Mixed-height radical example: J^(2) = (x^2) intersect (y,z)^2.")

# A nonmonomial, nonhomogeneous example with isolated primes of heights 1,3.
T, (a, b, c, d) = polynomial_ring(QQ, [:a, :b, :c, :d])
f = a^2 + b^2 - c
P = ideal(T, [f])
Q = ideal(T, [a - 1, b - c^2, d])
H = intersect(P, Q)
@assert associated_prime_heights(H) == [1, 3]
@assert symbolic_power(H, 2) == intersect(P^2, Q^2)
println("Nonhomogeneous mixed example: heights = ", associated_prime_heights(H))
