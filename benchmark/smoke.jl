using Oscar, SymbolicPowers

function timed(label, f)
  f() # Compile and warm caches before measuring.
  result = @timed f()
  println(label, ": ", result.time, " s, ", result.bytes, " bytes")
  return result.value
end

R, x = polynomial_ring(QQ, :x => 1:40)
triangle = ideal(R, [x[1] * x[2], x[1] * x[3], x[2] * x[3]])
timed("packing with 37 inactive variables", () -> is_packed(triangle))
timed("first packing witness", () -> non_packed_substitutions(triangle; first_only=true))
S, (a, b, c) = polynomial_ring(QQ, [:a, :b, :c])
I = ideal(S, [a * b, a * c, b * c])
timed("containment threshold, ordinary exponent 5", () -> containment_problem(I, 5))
timed("homogeneous initial degree", () -> initial_degree(I))
