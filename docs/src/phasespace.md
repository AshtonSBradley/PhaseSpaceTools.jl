# Usage Examples

Coherent states play a central role in quantum phase space methods, providing a mapping of many-body boson dynamics to equivalent stochastic differential equations.

## Glauber-Sudarshan P-representation
One of the most well-known representations is the Glauber-Sudarshan P distribution (here referred to as P). For a single mode optical field, the representation of the density matrix is given by

```math
\rho = \int d^2\alpha |\alpha\rangle\langle\alpha|P(\alpha,\alpha^*)
```
where the integral is taken over the entire complex plane. The P-representation occurs in a single complex variable phase space, as seen from the expansion over diagonal coherent state projectors.

In the P-representation a coherent state may be sampled as a single point on the complex plane

```julia
N = 10000
state = Coherent(12.0)
α,ᾱ = glauberP(state,N)
```

In any other distribution the coherent state is more interesting to sample.

## Positive-P
In the Positive-P representation the simples way to sample the state is again as a point on the complex plane

```julia
α,ᾱ = positiveP(state,N)
```

Here ``α,ᾱ`` are once again complex conjugates.
However, we emphasize that this is no longer the case in any subsequent dynamical evolution of the ensemble of phase space points. Furthermore, unlike the Glauber-P distribution, the positive-P distribution for a particular quantum state is not unique.

In this package we primarily focus on single mode problems, taking the general approach of sampling the Husimi-Q funcion, and then exploiting a simple convolutional relationship between the Husimi-Q and a particular choice of positive-P representation.

Once choice of positive-P that can be made for any single mode density matrix is given by
```math
P(\alpha,\bar\alpha) = \frac{1}{4\pi^2}\langle (\alpha+\bar\alpha^*)/2|\rho|(\alpha+\bar\alpha^*)/2\rangle e^{-|\alpha-\bar\alpha^*|^2/4}
```
We can then use the transformation
```math
\mu = (\alpha+\bar\alpha^*)/2,\quad\quad\gamma = (\alpha-\bar\alpha^*)/2,
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
is the [Husimi-Q function](https://en.wikipedia.org/wiki/Husimi_Q_representation), and the remaining complex Gaussian distribution is readily sampled as
```math
\gamma = \frac{1}{\sqrt{2}}\left(\eta_1+i\eta_2\right).
```
Here ``\eta_j`` are independent normal random variates with zero mean and unit variance:

```math
\langle \eta_j \rangle = 0\quad\quad \langle \eta_i\eta_j\rangle =\delta_{ij}.
```
Normally distributed Gaussian variates are sampled using `randn()`.
## Fock state
A more difficult state to sample an eigenstate of the number operator.
