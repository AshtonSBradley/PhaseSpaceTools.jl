#!/usr/bin/env julia

using PhaseSpaceTools
using Base.Test

#test coherent
include("samplecoherentW.jl")
@testset "coherent state W tests " begin include("testcoherentW.jl") end

# test fock
include("samplefockW.jl")
@testset "fock state W tests " begin include("testfockW.jl") end

include("samplefock+W.jl")
@testset "fock state +W tests " begin include("testfock+W.jl") end

include("samplefock+Wasymp.jl")
@testset "fock state +W tests for large n " begin include("testfock+Wasymp.jl") end

#test squeezed
include("samplesqueezed+P.jl")
@testset "squeezed state +P tests " begin include("testsqueezed+P.jl") end

#test squeezed
include("samplesqueezedW.jl")
@testset "squeezed state W tests " begin include("testsqueezed+P.jl") end


#test thermal
#include("samplethermalP.jl")
