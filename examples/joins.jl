# Elimination joins and their relation to symbolic powers.
using Oscar, SymbolicPowers

R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
A = ideal(R, [x, y])
B = ideal(R, [x, z])
@assert symbolic_join(A, B) == ideal(R, [x])
println("The join of the y-axis and z-axis is the plane x = 0.")

I = ideal(R, [x * y, x * z, y * z])
m = ideal(R, [x, y, z])
expected = I^2 + ideal(R, [x * y * z])
@assert symbolic_join(I, m^2) == expected
@assert symbolic_power_join(I, 2) == expected

# I is an intersection of the three coordinate primes, so radicality is
# independently known. Only in that situation is skipping its check justified.
@assert symbolic_power_join(I, 2; check=false) == expected
println("For the radical triangle ideal, I star m^2 = I^(2).")

# The general join is defined for nonradical ideals too, but is then a
# different operation: the symbolic-power helper checks this hypothesis.
N = ideal(R, [x^2])
@assert symbolic_join(N, m^2) == ideal(R, [x^3])
@assert symbolic_power(N, 2) == ideal(R, [x^4])
println("Nonradical warning: (x^2) star m^2 = (x^3), whereas (x^2)^(2) = (x^4).")
