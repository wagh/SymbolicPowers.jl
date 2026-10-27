using Oscar, SymbolicPowers
R, (x, y) = polynomial_ring(QQ, [:x, :y])
I = ideal(R, [x^3, x * y])
S = symbolic_power(I, 2; use_minimal_primes=true)
println("CHECK|minimal_square_is_x2|", S == ideal(R, [x^2]))
println("CHECK|actual_minimal_containment|", is_subset(S, I))
println(
  "CHECK|reported_minimal_containment|",
  is_symbolic_power_contained_in_power(I, 2, 1; use_minimal_primes=true),
)
println(
  "CHECK|reported_minimal_reverse_containment|",
  containment_problem(I, 2; in_symbolic=true, use_minimal_primes=true),
)
R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
I = ideal(R, [x * y, x * z, y * z])
println("CHECK|zeroth_power_is_unit|", is_one(symbolic_power(I, 0)))
try
  symbolic_power_join(ideal(R, [x^2]), 2)
  println("CHECK|nonradical_join_rejected|false")
catch err
  println("CHECK|nonradical_join_rejected|", err isa ArgumentError)
end
R, (x, y, z, w) = polynomial_ring(QQ, [:x, :y, :z, :w])
J = ideal(R, [(x + y) * y, (x + y) * z, y * z])
I = intersect(J, ideal(R, gens(R))^3)
println("CHECK|embedded_default_symbolic_square_is_ordinary|", symbolic_power(I, 2) == I^2)
println("CHECK|old_associated_heights|", associated_prime_heights(I))
println("CHECK|new_square_associated_heights|", associated_prime_heights(I^2))
println("CORRECTNESS_COMPLETE")
