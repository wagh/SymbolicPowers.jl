# Frobenius powers differ from ordinary powers; symbolic Frobenius helpers
# apply to primes. Small characteristic is also handled by Hasse derivatives.
using Oscar, SymbolicPowers

R, (x, y, z) = polynomial_ring(GF(5), [:x, :y, :z])
L = ideal(R, [x, y])
@assert frobenius_power(L, 5) == ideal(R, [x^5, y^5])
@assert x^4 * y in L^5
@assert !(x^4 * y in frobenius_power(L, 5))
println("In characteristic 5, (x,y)^[5] = (x^5,y^5) is strictly smaller than (x,y)^5.")

# This prime is a complete intersection, so its symbolic powers are ordinary.
Q = ideal(R, [x + y, z])
@assert symbolic_power_frobenius(Q, 3) == Q^3
@assert symbolic_power(Q, 3; method=:frobenius) == Q^3
println("Complete-intersection prime: the Frobenius helper returns Q^3 directly.")

# A non-complete-intersection prime: the kernel of x->t^3,y->t^4,z->t^5.
P = ideal(R, [y^2 - x * z, x^2 * y - z^2, x^3 - y * z])
S = symbolic_power_frobenius(P, 3)
@assert S == symbolic_power_join(P, 3)
println("Monomial curve over GF(5): Frobenius and join methods agree at order 3.")

# In characteristic 2, ordinary derivatives miss the order of a square.
# The join implementation uses Hasse derivatives, and works even for n > p.
T, (u, v) = polynomial_ring(GF(2), [:u, :v])
H = ideal(T, [u + v])
@assert symbolic_power_join(H, 3) == H^3
@assert symbolic_power_frobenius(H, 3) == H^3
@assert !((u + v)^2 in symbolic_power_join(H, 3))
println("Characteristic 2: (u+v)^2 is correctly excluded from the third symbolic power.")
