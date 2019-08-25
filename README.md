# PhaseSpaceTools

<!-- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://AshtonSBradley.github.io/PhaseSpaceTools.jl/stable) -->
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://AshtonSBradley.github.io/PhaseSpaceTools.jl/dev)
[![Build Status](https://travis-ci.com/AshtonSBradley/PhaseSpaceTools.jl.svg?branch=master)](https://travis-ci.com/AshtonSBradley/PhaseSpaceTools.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/t6i7kdnpffgnq6pg?svg=true)](https://ci.appveyor.com/project/AshtonSBradley/phasespacetools-jl)
[![Coverage Status](https://coveralls.io/repos/github/AshtonSBradley/PhaseSpaceTools.jl/badge.svg?branch=master)](https://coveralls.io/github/AshtonSBradley/PhaseSpaceTools.jl?branch=master)
[![codecov](https://codecov.io/gh/AshtonSBradley/PhaseSpaceTools.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/AshtonSBradley/PhaseSpaceTools.jl)

Package for sampling some of the quantum initial states commonly encountered in quantum-optical and matter-wave bosonic phase space simulations. Wigner and positive-P distributions are available, being the most useful for dynamical simulations.

Available distributions are `glauberP`, `positiveP` `wigner`, `positiveW`, `husimiQ`.

## Install

```julia
] add https://github.com/AshtonSBradley/PhaseSpaceTools.jl.git
```

To get help type, e.g.

```julia
julia> methods(positiveP)
```
for a list of available methods.

## Usage
```julia
julia> using PhaseSpaceTools
help?> positiveP

  search: positiveP positiveW

    α,ᾱ = positiveP(state <: State,N)

    Generate N samples from the positive-P phase-space distribution for state.

    Moments of the positive-P distribution generate quantum operator averages that are normally ordered.
```
## Implemented states

* Coherent(α)
* Thermal(α,n̄)
* Squeezed(α,ϵ)
* Fock(N)
* Bogoliubov(u,v,n̄)
* Crescent(α,ϵ,q)

### Coherent state
A coherent state |α⟩ is sampled as
```julia
α = 1.0+im*2.0 #coherent state amplitude
state = Coherent(α) # create state |α⟩
N = 1000 #number of samples
α,ᾱ = positiveP(state,N)
```
This is a special case where the two phase space variables `\alpha` and `\bar\alpha` are complex conjugate, and non-stochastic in the `+P` representation.

### Fock state
An approximate fock state sampler in the Wigner representation:
```julia
n = 100
state = Fock(n)  
N = 1000 #number of samples
a,ā = wigner(state,N)
```
Provides an approximate sampling of `W` that reproduces operator averages for large `n`.

## Examples

See  `/examples/PhaseSpaceTools.ipynb` for more usage.

# External links
___Numerical representation of quantum states in the positive-P and Wigner representations,___ \
M K Olsen, A S Bradley, \
[Optics Communications 282, 3924 (2009)](http://dx.doi.org/10.1016/j.optcom.2009.06.033)

[Scipost Commentary with full erratum](https://scipost.org/commentaries/10.1016/j.optcom.2009.06.033/)
