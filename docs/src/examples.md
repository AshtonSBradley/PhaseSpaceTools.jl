# Usage Examples

## Coherent state
A trivial example is provided by the coherent state $|\alpha\rangle$. These are eigenstate of the Bose annihilation operator ``a``, with commutation relation ``[a,a^\dagger]=1``. They satisfy the eigenvalue equation

```math
 a|\alpha\rangle = \alpha|\alpha\rangle
```

The Fock basis ``|n\rangle`` is often convenient to work in. Fock states are eigenstates of the number operator ``n=a^\dagger a``:

```math
a^\dagger a|n\rangle = n|n\rangle.
```

In the Fock basis, the coherent states take the form

```math
|\alpha\rangle = e^{-|\alpha|^2/2}\sum_{n=0}^\infty\frac{\alpha^n}{n!}|n\rangle.
```

Coherent states play a central role in quantum phase space methods, providing a mapping of many-body boson dynamical problems to equivalent stochastic differential equations.

### Glauber-P
In the Glauber-P representation the state may be sampled as a single point on the complex plane

```julia
N = 10000
state = Coherent(12.0)
α,ᾱ = glauberP(state,N)
```
### Positive-P
In the Positive-P representation the simples way to sample the state is again as a point on the complex plane

```julia
α,ᾱ = positiveP(state,N)
```

## Fock state
A more difficult state to sample is the eigenstate of the number operator...
