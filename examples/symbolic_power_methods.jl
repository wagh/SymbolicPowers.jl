# Compare all characteristic-zero methods on an ideal satisfying their hypotheses.
# The three coordinate points in P^2 have I^(2) = I^2 + (xyz).
using Oscar, SymbolicPowers

R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
I = ideal(R, [x * y, x * z, y * z])
expected = I^2 + ideal(R, [x * y * z])

# This radical ideal is monomial, unmixed and homogeneous of dimension one.
# Its primary components (x,y), (x,z), (y,z) are complete intersections.
available_methods = (:auto, :primary, :monomial, :saturation, :unmixed,
  :complete_intersection, :join)
for method in available_methods
  S = symbolic_power(I, 2; method)
  @assert S == expected
  println(method, ": I^(2) = I^2 + (xyz)")
end

@assert x * y * z in expected
@assert !(x * y * z in I^2)
@assert !is_symbolic_equal_ordinary(I, 2)
println("The element xyz witnesses the strict inclusion I^2 < I^(2).")

# In contrast, a complete-intersection prime has ordinary symbolic powers.
P = ideal(R, [x - y, y - z])
for n in 1:3
  @assert symbolic_power(P, n; method=:complete_intersection) == P^n
end
println("For the point ideal (x-y,y-z), symbolic and ordinary powers agree.")
