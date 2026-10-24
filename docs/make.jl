using Documenter
using Oscar
using SymbolicPowers

DocMeta.setdocmeta!(SymbolicPowers, :DocTestSetup,
  :(using Oscar, SymbolicPowers); recursive=true)

makedocs(;
  modules=[SymbolicPowers],
  sitename="SymbolicPowers.jl",
  format=Documenter.HTML(; prettyurls=get(ENV, "CI", "false") == "true",
    edit_link=nothing, repolink=nothing),
  remotes=nothing,
  checkdocs=:exports,
  pages=[
    "Introduction" => "index.md",
    "Mathematical guide" => "mathematics.md",
    "Algorithms and guarantees" => "algorithms.md",
    "API" => "api.md",
    "Development and validation" => "development.md",
  ],
)
