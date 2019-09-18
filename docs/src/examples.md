# Usage Examples

Coherent states play a central role in quantum phase space methods, providing a mapping of many-body boson dynamics to equivalent stochastic differential equations.

### Glauber-P
In the Glauber-P representation the state may be sampled as a single point on the complex plane

```julia
N = 10000
state = Coherent(12.0)
α,ᾱ = glauberP(state,N)
```

In any other distribution the coherent state is more interesting to sample.
### Positive-P
In the Positive-P representation the simples way to sample the state is again as a point on the complex plane

```julia
α,ᾱ = positiveP(state,N)
```

Here ``α,ᾱ`` are once again complex conjugates.
However, we emphasize that this is no longer the case in any subsequent dynamical evolution of the ensemble of phase space points. Furthermore, unlike the Glauber-P distribution, the positive-P distribution for a particular quantum state is not unique.

In this package, we treat single mode problems, taking the general approach of sampling the Husimi-Q funcion, and then exploiting the simple convolutional relationship between the Husimi-Q and the positive-P representations:

```math
P(\alpha,\bar\alpha) = \frac{1}{4\pi^2}\langle (\alpha+\bar\alpha^*)/2|\rho|(\alpha+\bar\alpha^*)/2\rangle e^{-|\alpha-\bar\alpha^*|^2/4}
```
We use the transformation
```math
\mu = (\alpha+\bar\alpha^*)/2,\quad\quad\gamma = (\alpha-\bar\alpha^*)/2
```
with inverse
```math
\alpha = \mu + \gamma,\quad\quad \bar\alpha = \mu^* - \gamma^*.
```
The problem reduces to that of sampling the  distribution
```math
P(\mu,\gamma)=Q(\mu,\mu^*)\frac{e^{-|\gamma|^2}}{\pi}
```
where
```math
Q(\alpha,\alpha^*)\equiv \frac{\langle \alpha|\rho|\alpha\rangle}{\pi}
```
is the Husimi-Q function, and the remaining Gaussian is readily sampled.
## Fock state
A more difficult state to sample an eigenstate of the number operator.
