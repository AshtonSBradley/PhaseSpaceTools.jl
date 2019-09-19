# Overview
This package provides sampling methods for commonly used quantum states in various quantum phase-space representations. These include Glauber-Sudarshan-P, Positive-P, Husimi-Q, and Wigner representations.

There are also convenience methods for calculating operator averages from phase-space averages, and for sampling noises for solving SDEs using [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl).

# Installation
In the julia REPL

```
julia> ]add https://github.com/AshtonSBradley/PhaseSpaceTools.jl

Pkg> test PhaseSpaceTools
```
# Reference
If you find this package useful, please cite

M. K. Olsen, A. S. Bradley, Numerical representation of quantum states in the positive-P and Wigner representations, [Optics Communications __282__, 3924-3929 (2009).](https://dx.doi.org/10.1016/j.optcom.2009.06.033)

A full erratum is given in this [SciPost Commentary.](https://scipost.org/commentaries/10.1016/j.optcom.2009.06.033/)
