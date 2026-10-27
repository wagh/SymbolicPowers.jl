# Symbolic defect counts missing generators; regularity data are finite samples.
using Oscar, SymbolicPowers

R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
I = ideal(R, [x * y, x * z, y * z])
@assert symbolic_power(I, 2) == I^2 + ideal(R, [x * y * z])
@assert symbolic_defect(I, 1) == 0
@assert symbolic_defect(I, 2) == 1
@assert min_degree_symbolic_power(I, 2) == 3
println("Triangle: defect at order 2 is 1, represented by the class of xyz.")

P = ideal(R, [x, y])
for n in 1:3
  @assert symbolic_defect(P, n) == 0
end
samples = symbolic_regularity_samples(P; sample_size=3)
@assert samples == [(power=n, regularity=n, ratio=QQ(1)) for n in 1:3]
@assert asymptotic_regularity(P; sample_size=3, approximate=true) == 1
println("Linear complete-intersection samples: ", samples)

# A less uniform family: the first three ratios are 4, 9/2, 13/3.
# Taking their minimum does not determine an asymptotic invariant or bound it.
T, (a, b, c, p, q, r, s, t, u) = polynomial_ring(QQ,
  [:a, :b, :c, :p, :q, :r, :s, :t, :u])
J = ideal(T, [a * b * c, b * c * p * q, a * c * r * s, a * b * t * u])
data = symbolic_regularity_samples(J; sample_size=3)
@assert [row.regularity for row in data] == [4, 9, 13]
statistic = asymptotic_regularity(J; sample_size=3, approximate=true)
@assert statistic == 4
println("Nonconstant regularity ratios: ", [row.ratio for row in data])
println("Their minimum, ", statistic, ", is only a finite statistic.")
