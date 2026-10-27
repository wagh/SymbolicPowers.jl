# Containment predicates and searches in both directions.
using Oscar, SymbolicPowers

R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
I = ideal(R, [x * y, x * z, y * z])
@assert !is_symbolic_power_contained_in_power(I, 2, 2)
@assert is_symbolic_power_contained_in_power(I, 3, 2)
@assert containment_problem(I, 2) == 3
@assert containment_problem(I, 3; in_symbolic=true) == 2
@assert !is_symbolic_equal_ordinary(I, 2)

println("For each n, find the least m with I^(m) contained in I^n:")
for n in 1:4
  m = containment_problem(I, n)
  @assert is_symbolic_power_contained_in_power(I, m, n)
  @assert !is_symbolic_power_contained_in_power(I, m - 1, n)
  println("  n = ", n, ", least m = ", m)
end
println("Reverse search: the largest k with I^(3) contained in I^k is 2.")

# Minimal-prime powers need a different search interval. The only minimal
# prime of J is (x), so J^(m)_min = (x^m), even when m = 1.
J = ideal(R, [x^3, x * y])
@assert symbolic_power(J, 2; use_minimal_primes=true) == ideal(R, [x^2])
@assert !is_symbolic_power_contained_in_power(J, 2, 1; use_minimal_primes=true)
@assert containment_problem(J, 2; in_symbolic=true, use_minimal_primes=true) == 0
@assert containment_problem(J, 2; use_minimal_primes=true, max_power=6) == 6
println("Minimal-prime example: (x^2) is not contained in J; the reverse answer is 0.")
println("The least m with (x^m) contained in J^2 is 6 (max_power bounds the search).")
