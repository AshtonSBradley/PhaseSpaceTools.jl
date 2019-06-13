#test
using Pkg
Pkg.activate(".")

using Test, LinearAlgebra, Revise, PhaseSpaceTools

Pkg.test("PhaseSpaceTools")
