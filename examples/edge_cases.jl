# Zero/unit ideals, exponent zero, zero-variable rings and rejected hypotheses.
using Oscar, SymbolicPowers

R, (x, y) = polynomial_ring(QQ, [:x, :y])
Z = ideal(R, elem_type(R)[])
U = ideal(R, [one(R)])
I = ideal(R, [x * y])
@assert is_one(symbolic_power(Z, 0))
@assert is_one(symbolic_power(I, 0))
for n in 1:3
  @assert is_zero(symbolic_power(Z, n))
  @assert is_one(symbolic_power(U, n))
end
@assert big_height(Z) == 0
@assert initial_degree(U) == 0
println("Zeroth symbolic powers are the unit ideal, including for the zero ideal.")

# A polynomial ring in no variables is just its coefficient field.
T, _ = polynomial_ring(QQ, Symbol[])
Z0 = ideal(T, elem_type(T)[])
U0 = ideal(T, [one(T)])
for A in (Z0, U0), B in (Z0, U0)
  S = symbolic_join(A, B)
  @assert is_zero(S) == (is_zero(A) && is_zero(B))
  @assert is_one(S) == (is_one(A) || is_one(B))
end
println("In a zero-variable polynomial ring, the join is the sum of the ideals.")

# Expected errors are caught explicitly so this file remains runnable.
# Failed hypotheses must not silently select an unrelated algorithm.
rejected = (
  ("negative exponent", () -> symbolic_power(I, -1)),
  ("unknown method", () -> symbolic_power(I, 2; method=:unknown)),
  ("nonradical join input", () -> symbolic_power_join(ideal(R, [x^2]), 2)),
  ("Frobenius in characteristic zero", () -> frobenius_power(I, 2)),
  ("initial degree of the zero ideal", () -> initial_degree(Z)),
  ("exact Waldschmidt constant of nonmonomial input",
    () -> waldschmidt_constant(ideal(R, [x - y]))),
  ("regularity statistic without explicit opt-in",
    () -> asymptotic_regularity(I; sample_size=2)),
)
for (description, compute) in rejected
  failure = try
    compute()
    nothing
  catch err
    err
  end
  @assert failure isa ArgumentError
  println("Rejected ", description, ": ", sprint(showerror, failure))
end
