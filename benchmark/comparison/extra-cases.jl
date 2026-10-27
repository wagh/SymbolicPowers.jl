const RF5, _ = polynomial_ring(GF(5), :k => 0:2)
const RF101, _ = polynomial_ring(GF(101), :l => 0:2)
const POINT_COORDINATES = [
  [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1], [1, 2, 3], [3, 1, -2], [2, 3, 1], [1, -1, 2]
]
function point_ideal(R, p)
  j = findfirst(!iszero, p)
  return ideal(R, [p[j] * R[i] - p[i] * R[j] for i in 1:3 if i != j])
end
eight_points() =
  ideal(R3, gens(intersect([point_ideal(R3, p) for p in POINT_COORDINATES]...)))
fat_points() = ideal(
  R3,
  gens(intersect([point_ideal(R3, POINT_COORDINATES[i])^(i == 1 ? 2 : 1) for i in 1:4]...)),
)
mixed_ideal() = ideal(
  R4,
  gens(
    intersect(
      ideal(R4, [R4[1]^2 + R4[2]^2 - R4[3]]),
      ideal(R4, [R4[1] - 1, R4[2] - R4[3]^2, R4[4]]),
    ),
  ),
)
const EXTRA_CASES = [
  ("points8_sp2", eight_points, I -> symbolic_power(I, 2), :ideal),
  ("points8_sp4", eight_points, I -> symbolic_power(I, 4), :ideal),
  ("fatpoints4_sp2", fat_points, I -> symbolic_power(I, 2), :ideal),
  (
    "points8_reg3",
    eight_points,
    I -> asymptotic_regularity(I; sample_size=3, approximate=true),
    :value,
  ),
  (
    "fermat_ci3",
    () -> fermat(R3),
    I -> symbolic_power(I, 3; method=:complete_intersection),
    :ideal,
  ),
  ("mixed_sp2", mixed_ideal, I -> symbolic_power(I, 2), :ideal),
  (
    "frobenius_linear101_sp2",
    () -> ideal(RF101, [RF101[1] + RF101[2], RF101[3]]),
    I -> symbolic_power_frobenius(I, 2),
    :ideal,
  ),
  (
    "frobenius_point5_sp2",
    () -> ideal(RF5, [RF5[1] - 1, RF5[2], RF5[3] + RF5[1]]),
    I -> symbolic_power_frobenius(I, 2),
    :ideal,
  ),
]
filter!(
  spec -> !(spec[1] in ("contain_triangle24", "squarefree_m7", "squarefree_m10")), CASES
)
append!(CASES, EXTRA_CASES)
const MONO_IDS = ("triangle_sp6", "triangle_sp12", "embedded_sp3", "ci_monomial_sp8",
  "contain_triangle8", "equal_triangle4", "wald_triangle", "wald_cycle5",
  "defect_triangle2",
  "resurgence_triangle5", "reg_triangle4", "konig_complete6", "konig_hyper6",
  "packed_cycle4", "packed_cycle5", "packed_triangle_ambient8", "squarefree_m5",
  "squarefree_m6")
append!(
  CASES,
  [
    ("mono_" * id, make, operation, kind) for
    (id, make, operation, kind) in CASES if id in MONO_IDS
  ],
)
