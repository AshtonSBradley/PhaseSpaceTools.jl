# PhaseSpaceTools

<!-- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://AshtonSBradley.github.io/PhaseSpaceTools.jl/stable) -->
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://AshtonSBradley.github.io/PhaseSpaceTools.jl/dev)
[![Build Status](https://travis-ci.com/AshtonSBradley/PhaseSpaceTools.jl.svg?branch=master)](https://travis-ci.com/AshtonSBradley/PhaseSpaceTools.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/t6i7kdnpffgnq6pg?svg=true)](https://ci.appveyor.com/project/AshtonSBradley/phasespacetools-jl)
[![Coverage Status](https://coveralls.io/repos/github/AshtonSBradley/PhaseSpaceTools.jl/badge.svg?branch=master)](https://coveralls.io/github/AshtonSBradley/PhaseSpaceTools.jl?branch=master)
[![codecov](https://codecov.io/gh/AshtonSBradley/PhaseSpaceTools.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/AshtonSBradley/PhaseSpaceTools.jl) [![DOI](https://zenodo.org/badge/115932136.svg)](https://zenodo.org/badge/latestdoi/115932136)



Package for sampling some of the quantum initial states commonly encountered in quantum-optical and matter-wave bosonic phase space simulations. Wigner and positive-P distributions are available, being the most useful for dynamical simulations.

Available distributions are `glauberP`, `positiveP` `wigner`, `positiveW`, `husimiQ`.

## Install

```julia
] add https://github.com/AshtonSBradley/PhaseSpaceTools.jl.git
```

## Usage
```julia-repl
julia> using PhaseSpaceTools

help?> positiveP
search: positiveP positiveW

  α,α⁺ = positiveP(state <: State,N)

  Generate N samples from the positive-P (+P) phase-space distribution for state.

  Moments of the +P distribution generate quantum operator averages that are normally ordered.

  In general the two random variates α,α⁺ are statistically independent for the +P distribution. 
```
## Implemented states

```julia-repl
help?> State
search: State state estimate InvalidStateException AbstractSet AbstractVector AbstractVecOrMat stacktrace StackTraces istaskstarted abstract type AbstractRange AbstractPattern

  State

  Abstract container type for all sampled states: state <: State.

  Examples
  ≡≡≡≡≡≡≡≡≡≡

  Find all states that may be sampled
  =====================================

  julia> subtypes(State)
  
  6-element Vector{Any}:
  Bogoliubov
  Coherent
  Crescent
  Fock
  Squeezed
  Thermal

  Create and sample a particular state (vacuum)
  ===============================================

  julia> s = Fock(0)
  Fock(0)
  julia> wigner(s,100)
  julia> wigner(s,100)
  ┌ Warning: Fock state sampling for W is only valid for n ≫ 1.

  Here a warning is generated since Fock sampling is not well defined for small n. 

  A simpler way to sample the vacuum is 

  julia> s = Coherent(0) 
  Coherent(0.0 + 0.0im)  # type conversion to ComplexF64.
  
  julia> wigner(s,100)
  (ComplexF64[0.33820868828162637 + 0.4407579103538181im, 0.057183146091823775 - 0.2772571883006981im, ...

  generating two vectors of sampled points α,α⁺ in the complex plane. In this case, α = conj(α⁺), as we are not working with a doubled phase space.
```

### Coherent state
A coherent state |α⟩ is sampled as
```julia-repl
α = 1.0+im*2.0 # coherent amplitude
s = Coherent(α) # define state |α⟩
N = 1000 # number of samples
a,a⁺ = positiveP(s,N)
```
This is a special case where the two phase space variables `a` and `a⁺` are complex conjugate, and non-stochastic in the `+P` representation.

### Fock state
An approximate Fock state sampler in the Wigner representation:
```julia-repl
n = 100
s = Fock(n) # define number state |n⟩
N = 1000 # number of samples
a,a⁺ = wigner(s,N)
```
Provides an approximate sampling of `W` that reproduces operator averages for large `n`.

## Examples

See  `/examples/PhaseSpaceTools.ipynb` for more usage.

# External links
___Numerical representation of quantum states in the positive-P and Wigner representations,___ \
M K Olsen, A S Bradley, \
[Optics Communications 282, 3924 (2009)](http://dx.doi.org/10.1016/j.optcom.2009.06.033)

[Scipost Commentary with full erratum](https://scipost.org/commentaries/10.1016/j.optcom.2009.06.033/)
