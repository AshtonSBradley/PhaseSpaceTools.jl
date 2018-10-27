#test juila 1.0
using Pkg
pkg"activate ."
using PhaseSpaceTools, Revise
using Test
pkg"test PhaseSpaceTools"
