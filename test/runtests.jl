#!/usr/bin/env julia

using PhaseSpaceTools
using Test

@testset "coherent state W" begin include("testcoherentW.jl") end

@testset "fock state W" begin include("testfockW.jl") end

@testset "fock state +W" begin include("testfock+W.jl") end

@testset "fock state +W for large n " begin include("testfock+Wasymp.jl") end

@testset "squeezed state +P " begin include("testsqueezed+P.jl") end

@testset "squeezed state W" begin include("testsqueezed+P.jl") end

@testset "Bogoliubov W" begin include("testBogoliubov.jl") end

#TODO
# @testset "Glauber P" begin include("testthermalP.jl") end
