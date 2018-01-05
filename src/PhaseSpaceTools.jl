__precompile__()

module PhaseSpaceTools

using PyCall, SymPy, Distributions

#import laguerre from scipy.special
const scipy_spec = PyNULL()

function __init__()
    copy!(scipy_spec, pyimport_conda("scipy.special", "scipy"))
end

#using Reexport
#@reexport using Distributions

include("reject.jl")
include("laguerre.jl")
include("plaguerre.jl")
include("coherent.jl")
include("thermal.jl")
include("squeezed.jl")
include("fock.jl")
include("crescent.jl")

export coherent, thermal, squeezed, fock, crescent

end #end of module
