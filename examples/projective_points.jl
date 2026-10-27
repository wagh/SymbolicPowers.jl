# Reduced points, fat points and an embedded irrelevant component in P^2.
using Oscar, SymbolicPowers

R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
P1 = ideal(R, [y, z])        # (1:0:0)
P2 = ideal(R, [x, z])        # (0:1:0)
P3 = ideal(R, [x, y])        # (0:0:1)
P4 = ideal(R, [x - z, y - z]) # (1:1:1)
I = intersect(P1, P2, P3, P4)

# At each point the symbolic square prescribes vanishing to order at least 2.
expected = intersect(P1^2, P2^2, P3^2, P4^2)
S = symbolic_power(I, 2; method=:saturation)
@assert S == expected
@assert symbolic_power(I, 2) == S
println(
  "Four reduced points: the symbolic square is the intersection of their squared ideals."
)

# Doubling the first point changes its order in the symbolic square to 4.
F = intersect(P1^2, P2, P3, P4)
fat_square = symbolic_power(F, 2)
@assert fat_square == intersect(P1^4, P2^2, P3^2, P4^2)
println("Fat points: multiplicities (2,1,1,1) become (4,2,2,2).")

# Adding an embedded component at m changes the default convention. The
# minimal-prime convention discards it and recovers the reduced-points answer.
m = ideal(R, [x, y, z])
E = intersect(I, m^3)
@assert associated_prime_heights(E) == [2, 3]
@assert symbolic_power(E, 2) == E^2
@assert symbolic_power(E, 2; use_minimal_primes=true) == S
println(
  "With an embedded irrelevant component: default E^(2) = E^2; minimal-prime E^(2) = I^(2)."
)
