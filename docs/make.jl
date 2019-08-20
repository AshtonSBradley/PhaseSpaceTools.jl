import Pkg; Pkg.add("Documenter")

using Documenter, PhaseSpaceTools

makedocs(sitename="PhaseSpaceTools.jl")

deploydocs(
    repo = "github.com/AshtonSBradley/PhaseSpaceTools.jl.git",
)
