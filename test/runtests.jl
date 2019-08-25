#!/usr/bin/env julia

using PhaseSpaceTools, Test
@testset "Testing methods" begin include("testmethods.jl") end

@testset "Coherent state W" begin include("testcoherentW.jl") end

@testset "Coherent state +P" begin include("testcoherent+P.jl") end

@testset "Fock state W" begin include("testfockW.jl") end

@testset "Fock state +W" begin include("testfock+W.jl") end

@testset "Fock state +W for large n " begin include("testfock+Wasymp.jl") end

@testset "Squeezed state +P " begin include("testsqueezed+P.jl") end

@testset "Squeezed state W" begin include("testsqueezed+P.jl") end

@testset "Bogoliubov W" begin include("testBogoliubov.jl") end

@testset "Noises" begin include("testNoises.jl") end

#TODO
# @testset "Glauber P" begin include("testthermalP.jl") end
