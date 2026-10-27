# Run from the package root: julia --project=. examples/run_all.jl
# Each file has its own module, so its variables and definitions are isolated.
let files = [
    "quickstart.jl",
    "symbolic_power_methods.jl",
    "embedded_and_mixed_ideals.jl",
    "projective_points.jl",
    "curves.jl",
    "joins.jl",
    "positive_characteristic.jl",
    "containment.jl",
    "polyhedra_and_waldschmidt.jl",
    "defect_and_regularity.jl",
    "packing.jl",
    "fields_and_gradings.jl",
    "edge_cases.jl",
  ]
  for file in files
    println("\n--- ", file, " ---")
    scope = Module(gensym(:SymbolicPowersExample))
    Base.include(scope, joinpath(@__DIR__, file))
  end
  println("\nAll ", length(files), " example files completed successfully.")
end
