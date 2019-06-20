# PhaseSpaceTools

[![Build Status](https://travis-ci.org/AshtonSBradley/PhaseSpaceTools.jl.svg?branch=master)](https://travis-ci.org/AshtonSBradley/PhaseSpaceTools.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/t6i7kdnpffgnq6pg?svg=true)](https://ci.appveyor.com/project/AshtonSBradley/phasespacetools-jl)
[![Coverage Status](https://coveralls.io/repos/AshtonSBradley/PhaseSpaceTools.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/AshtonSBradley/PhaseSpaceTools.jl?branch=master)
[![codecov.io](http://codecov.io/github/AshtonSBradley/PhaseSpaceTools.jl/coverage.svg?branch=master)](http://codecov.io/github/AshtonSBradley/PhaseSpaceTools.jl?branch=master)

Package for sampling some of the quantum initial states commonly encountered in quantum-optical and matter-wave bosonic phase space simulations. Wigner (`W`) and positive-P (`+P`) representations are implemented, being the most useful for dynamical simulations. Currently supports only single mode sampling.

Note: this package is being refactored at present, but should also now be fairly self documenting. 

Available distributions are `glauberP`, `wigner`, `husimiQ`, `positiveP`, `positiveW`

To get help type, e.g.

```julia
> methods(positiveP)

```

for a list of what is implemented

## Install

```julia
] add https://github.com/AshtonSBradley/PhaseSpaceTools.jl.git
```

## Usage
```julia
julia>using PhaseSpaceTools
help?> squeezed
search: squeezed

  a,ā = squeezed(β,ϵ,N;dist=:posP)

  Sample the phase-space distribution for a squeezed state.

  β: coherent (complex) amplitude.

  ϵ: complex valued squeezing parameter.

  N: number of samples.

  dist: phase-space distribution; can be :W or :posP.

  For standard P,Q,W distributions, a and ā are complex conjugate, while for +P etc, a and ā are independent
  variables.
```


### States
`coherent`, `thermal`, `squeezed`, `fock`, `crescent`

#### Coherent state
A coherent state |α⟩ is "sampled" as
```julia
α = 1.0+im*2.0 #coherent state amplitude
state = Coherent(α) # create state |α⟩
N = 1000 #number of samples
a,ā = positiveP(state,N)
```
This is a special case where the two phase space variables `a` and `ā` are complex conjugate, and non-stochastic in the `+P` representation.

#### Fock state
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
