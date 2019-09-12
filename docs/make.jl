import Pkg; Pkg.add("Documenter")

using Documenter, PhaseSpaceTools

makedocs(sitename="PhaseSpaceTools.jl",
pages = [
"Examples" => "examples.md",
"States" => "states.md",
"Distributions" => "distributions.md",
"Solving SDES" => "sdes.md"
])

deploydocs(
    repo = "github.com/AshtonSBradley/PhaseSpaceTools.jl.git",
)
