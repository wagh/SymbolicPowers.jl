# SymbolicPowers.jl: identical inputs to bench.m2; warm code, fresh ideals.
using Oscar, SymbolicPowers

println("META|Julia|", VERSION)
println("META|Oscar|", pkgversion(Oscar))
println("META|SymbolicPowers|", pkgversion(SymbolicPowers))

const R3, _ = polynomial_ring(QQ, [:x, :y, :z])
const RW, _ = graded_polynomial_ring(QQ, [:u, :v, :w]; weights=[3, 4, 5])
const R4, _ = polynomial_ring(QQ, [:a, :b, :c, :d])
const R5, _ = polynomial_ring(QQ, :e => 0:4)
const R6, _ = polynomial_ring(QQ, :f => 0:5)
const R7, _ = polynomial_ring(QQ, :g => 0:6)
const R8, _ = polynomial_ring(QQ, :h => 0:7)
const R10, _ = polynomial_ring(QQ, :j => 0:9)
const RF, _ = polynomial_ring(GF(3), [:p, :q, :r])

triangle(R) = ideal(R, [R[1] * R[2], R[1] * R[3], R[2] * R[3]])
function fermat(R)
  x, y, z = gens(R)
  return ideal(R, [x * (y^3 - z^3), y * (z^3 - x^3), z * (x^3 - y^3)])
end
function curve(R)
  x, y, z = gens(R)
  return ideal(R, [y^2 - x * z, x^2 * y - z^2, x^3 - y * z])
end
twisted(R) =
  ideal(R, [R[1] * R[3] - R[2]^2, R[1] * R[4] - R[2] * R[3], R[2] * R[4] - R[3]^2])
embedded(R) = ideal(R, [R[1]^3, R[1]^2 * R[2]^2, R[1] * R[3]^3, R[2]^4, R[2]^2 * R[3]])
cycle(R) = ideal(R, [R[i] * R[mod1(i + 1, ngens(R))] for i in 1:ngens(R)])
hyper(R) = ideal(
  R,
  [
    R[i] * R[j] * R[k] for i in 1:ngens(R) for j in (i + 1):ngens(R) for
    k in (j + 1):ngens(R)
  ],
)

const CASES = [
  ("triangle_sp6", () -> triangle(R3), I -> symbolic_power(I, 6), :ideal),
  ("triangle_sp12", () -> triangle(R3), I -> symbolic_power(I, 12), :ideal),
  ("fermat_sp2", () -> fermat(R3), I -> symbolic_power(I, 2), :ideal),
  ("fermat_sp3", () -> fermat(R3), I -> symbolic_power(I, 3), :ideal),
  ("curve_sp3", () -> curve(R3), I -> symbolic_power(I, 3), :ideal),
  ("curve_weighted_sp3", () -> curve(RW), I -> symbolic_power(I, 3), :ideal),
  ("twisted_sp3", () -> twisted(R4), I -> symbolic_power(I, 3), :ideal),
  ("embedded_sp3", () -> embedded(R3), I -> symbolic_power(I, 3), :ideal),
  (
    "embedded_min_sp3",
    () -> embedded(R3),
    I -> symbolic_power(I, 3; use_minimal_primes=true),
    :ideal,
  ),
  (
    "ci_monomial_sp8",
    () -> ideal(R3, [R3[1]^2, R3[2]^2]),
    I -> symbolic_power(I, 8),
    :ideal,
  ),
  (
    "affine_sp3",
    () -> ideal(R4, [R4[2] * R4[3] + 1, R4[1] * R4[3] * R4[4]]),
    I -> symbolic_power(I, 3),
    :ideal,
  ),
  ("join_triangle4", () -> triangle(R3), I -> symbolic_power_join(I, 4), :ideal),
  (
    "join_triangle4_core",
    () -> triangle(R3),
    I -> symbolic_power_join(I, 4; check=false),
    :ideal,
  ),
  (
    "join_twisted3_core",
    () -> twisted(R4),
    I -> symbolic_power_join(I, 3; check=false),
    :ideal,
  ),
  ("frobenius_curve3", () -> curve(RF), I -> symbolic_power_frobenius(I, 3), :ideal),
  ("frobenius_curve4", () -> curve(RF), I -> symbolic_power_frobenius(I, 4), :ideal),
  ("contain_triangle8", () -> triangle(R3), I -> containment_problem(I, 8), :value),
  ("contain_triangle24", () -> triangle(R3), I -> containment_problem(I, 24), :value),
  ("contain_fermat2", () -> fermat(R3), I -> containment_problem(I, 2), :value),
  ("equal_triangle4", () -> triangle(R3), I -> is_symbolic_equal_ordinary(I, 4), :value),
  ("height_fermat3", () -> fermat(R3), I -> big_height(I), :value),
  ("wald_triangle", () -> triangle(R3), I -> waldschmidt_constant(I), :value),
  ("wald_cycle5", () -> cycle(R5), I -> waldschmidt_constant(I), :value),
  ("defect_triangle2", () -> triangle(R3), I -> symbolic_defect(I, 2), :value),
  (
    "resurgence_triangle5",
    () -> triangle(R3),
    I -> lower_bound_resurgence(I; sample_size=5, use_waldschmidt=true),
    :value,
  ),
  (
    "reg_triangle4",
    () -> triangle(R3),
    I -> asymptotic_regularity(I; sample_size=4, approximate=true),
    :value,
  ),
  ("konig_complete6", () -> ideal(R6, gens(R6)), I -> is_konig(I), :value),
  ("konig_hyper6", () -> hyper(R6), I -> is_konig(I), :value),
  ("packed_cycle4", () -> cycle(R4), I -> is_packed(I), :value),
  ("packed_cycle5", () -> cycle(R5), I -> is_packed(I), :value),
  ("packed_triangle_ambient8", () -> triangle(R8), I -> is_packed(I), :value),
  ("squarefree_m5", () -> ideal(R5, gens(R5)), I -> squarefree_in_codim(I), :monomials),
  ("squarefree_m6", () -> ideal(R6, gens(R6)), I -> squarefree_in_codim(I), :monomials),
  ("squarefree_m7", () -> ideal(R7, gens(R7)), I -> squarefree_in_codim(I), :monomials),
  ("squarefree_m10", () -> ideal(R10, gens(R10)), I -> squarefree_in_codim(I), :monomials),
]

# Macaulay2 outputs generators in a coefficient/exponent interchange format.
# Compare actual ideals; equal generator order, GB order, or normalization is
# not assumed. This comparison is outside every timed region.
function parse_rational(s)
  a = split(strip(replace(s, "//" => "/")), '/')
  return QQ(parse(BigInt, a[1]), length(a) == 1 ? BigInt(1) : parse(BigInt, a[2]))
end
function read_m2_results(path)
  polys = Dict{String,Vector{Vector{Tuple{QQFieldElem,Vector{Int}}}}}()
  values = Dict{String,String}()
  for line in eachline(path)
    a = split(line, '|')
    if a[1] == "IDEAL"
      polys[a[2]] = Vector{Tuple{QQFieldElem,Vector{Int}}}[]
    elseif a[1] == "POLY"
      push!(polys[a[2]], Tuple{QQFieldElem,Vector{Int}}[])
    elseif a[1] == "TERM"
      exps = parse.(Int, split(replace(a[5], r"[{}\[\] ]" => ""), ','))
      push!(polys[a[2]][parse(Int, a[3]) + 1], (parse_rational(a[4]), exps))
    elseif a[1] == "VALUE"
      values[a[2]] = a[3]
    end
  end
  return polys, values
end

isempty(ARGS) && error("Pass the matching Macaulay2 output log as the first argument")
const M2_POLYS, M2_VALUES = read_m2_results(ARGS[1])

function agrees_with_m2(id, input, answer, kind)
  if kind == :value
    haskey(M2_VALUES, id) || return "missing"
    s = M2_VALUES[id]
    expected = if s == "true"
      true
    elseif s == "false"
      false
    else
      parse_rational(s)
    end
    return answer == expected ? "true" : "false"
  end
  haskey(M2_POLYS, id) || return "missing"
  R = base_ring(input)
  K = coefficient_ring(R)
  x = gens(R)
  gs = elem_type(R)[]
  for terms in M2_POLYS[id]
    f = zero(R)
    for (c, e) in terms
      coeff = K(numerator(c)) / K(denominator(c))
      f += coeff * prod(x[i]^e[i] for i in eachindex(x))
    end
    push!(gs, f)
  end
  expected = ideal(R, gs)
  actual = kind == :monomials ? ideal(R, answer) : answer
  return actual == expected ? "true" : "false"
end

function run_case(id, make_ideal::B, operation::F, kind) where {B,F}
  println("START|", id)
  flush(stdout)
  operation(make_ideal()) # Compile this operation before taking measurements.
  samples = Float64[]
  answer = nothing
  input = make_ideal()
  for _ in 1:5
    input = make_ideal()
    t = time_ns()
    answer = operation(input)
    push!(samples, (time_ns() - t) / 1e9)
  end
  ordered = sort(samples)
  println(
    "BENCH|", id, "|", ordered[3], "|", first(ordered), "|", last(ordered), "|", samples
  )
  matches = agrees_with_m2(id, input, answer, kind)
  println("MATCH|", id, "|", matches)
  @assert matches == "true" "Macaulay2 result is missing or differs for $id"
  flush(stdout)
end

include(joinpath(@__DIR__, "extra-cases.jl"))

if get(ENV, "BENCH_DEFINITIONS_ONLY", "0") != "1"
  selected = get(ENV, "BENCH_CASE", "")
  failures = String[]
  for spec in CASES
    isempty(selected) || spec[1] == selected || continue
    try
      run_case(spec...)
    catch err
      push!(failures, spec[1])
      println("ERROR|", spec[1], "|", sprint(showerror, err))
      flush(stdout)
    end
  end
  isempty(failures) || error("Benchmark failures: $(join(failures, ", "))")
  println("BENCHMARKS_COMPLETE")
  include(joinpath(@__DIR__, "correctness.jl"))
end
