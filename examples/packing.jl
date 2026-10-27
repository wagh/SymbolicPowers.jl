# Graph edge ideals: squarefree generators, disjoint supports and packing minors.
using Oscar, SymbolicPowers

R, v = polynomial_ring(QQ, :v => 1:6)
cycle4 = ideal(R, [v[1] * v[2], v[2] * v[3], v[3] * v[4], v[4] * v[1]])
@assert codim(cycle4) == 2
@assert is_konig(cycle4)
@assert is_packed(cycle4)
@assert squarefree_in_codim(cycle4) == [prod(v[1:4])]
println("The four-cycle is packed; its squarefree product in codimension is ",
  squarefree_in_codim(cycle4), ".")

triangle = ideal(R, [v[1] * v[2], v[1] * v[3], v[2] * v[3]])
@assert !is_konig(triangle)
@assert !is_packed(triangle)
witnesses = non_packed_substitutions(triangle)
@assert length(witnesses) == 1
@assert isempty(witnesses[1].zero) && isempty(witnesses[1].one)
println("The triangle itself is a non-Konig witness: ", witnesses[1])

# Attach a leaf to each triangle vertex. The three leaf edges are disjoint,
# so this graph is Konig, but deleting all leaves leaves a non-Konig triangle.
J = ideal(R, [v[1] * v[2], v[1] * v[3], v[2] * v[3],
  v[1] * v[4], v[2] * v[5], v[3] * v[6]])
@assert is_konig(J)
@assert !is_packed(J)
all_witnesses = non_packed_substitutions(J)
@assert any(w -> Set(w.zero) == Set(v[4:6]) && isempty(w.one), all_witnesses)
w = only(non_packed_substitutions(J; first_only=true))
println("A Konig ideal can fail packing; first witness: ", w)

# Verify a returned witness using an ordinary polynomial-ring map.
coordinate_images = [
  if a in w.zero
    zero(R)
  elseif a in w.one
    one(R)
  else
    a
  end for a in v
]
substitution_map = hom(R, R, coordinate_images)
minor_ideal = ideal(R, [substitution_map(g) for g in gens(J)])
@assert !is_konig(minor_ideal)

# The Konig predicate permits nonsquarefree monomials; packing requires a
# squarefree ideal. squarefree_gens filters the minimal monomial generators.
N = ideal(R, [v[1]^2 * v[2], v[1] * v[3], v[4]^3])
@assert squarefree_gens(N) == [v[1] * v[3]]
@assert is_konig(N)
println("Squarefree generators of a nonsquarefree ideal: ", squarefree_gens(N))
