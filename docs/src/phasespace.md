# Phase Space
For systems of Bose particles, the most important states are [coherent states](https://en.wikipedia.org/wiki/Coherent_states). These states allow representation of a large class of dynamical quantum evolution in terms of equivalent stochastic processes.

## Fock states
Coherent states may be represented in terms of [Fock states](https://en.wikipedia.org/wiki/Fock_state). For a single mode, these are just eigenstates of the number operator $a^\dagger a$, satisfying
```math
a^\dagger a|n\rangle = n|n\rangle
```
and these states are orthonormal $\langle n|m\rangle = \delta_{nm}$.

## Coherent states
[Coherent states](https://en.wikipedia.org/wiki/Coherent_states) play a central role in quantum phase space methods, providing a mapping of many-body boson dynamics to equivalent stochastic differential equations.

For a single bosonic mode, with creation and destruction operators $a^\dagger$, and $a$, with commutator $[a,a^\dagger]=1$, the coherent states are eigenstates of the annihilation operator
```math
a|\alpha\rangle =\alpha|\alpha\rangle.
```
where $\alpha\in \mathbb{C}$ can vary over the complex plane and plays the role of a classical amplitude.

In terms of Fock states, the coherent states take the form
```math
|\alpha\rangle = e^{-|\alpha|^2/2}\sum_{n=0}^\infty \frac{\alpha^n}{\sqrt{n!}}|n\rangle
```
These states are overcomplete, with resolution of the identity
```math
\mathbb{1} = \frac{1}{\pi}\int d^2\alpha |\alpha\rangle\langle\alpha|.
```
This property is responsible for much of the practical utility of coherent states, furnishing a mapping from quantum operators in Hilbert space to differential operators in phase space.

## Glauber-Sudarshan P-representation
One of the most well-known representations is the [Glauber-Sudarshan P-representation](https://en.wikipedia.org/wiki/Glauber–Sudarshan_P_representation) (for brevity here called the Glauber-P). For a single mode optical field, the representation of the density matrix is given by

```math
\rho = \int d^2\alpha |\alpha\rangle\langle\alpha|P(\alpha,\alpha^*)
```
where the integral is taken over the entire complex plane. The Glauber-P representation occurs in a single complex variable phase space, as seen from the expansion over diagonal coherent state projectors.

In the P-representation a coherent state may be sampled as a single point on the complex plane

```julia
N = 10000
state = Coherent(12.0)
α,ᾱ = glauberP(state,N)
```

In any other distribution the coherent state is more interesting to sample.

## Positive-P representation
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
is the [Husimi-Q representation](https://en.wikipedia.org/wiki/Husimi_Q_representation), and the remaining complex Gaussian distribution is readily sampled as
```math
\gamma = \frac{1}{\sqrt{2}}\left(\eta_1+i\eta_2\right).
```
Here ``\eta_j`` are independent normal random variates with zero mean and unit variance:

```math
\langle \eta_j \rangle = 0\quad\quad \langle \eta_i\eta_j\rangle =\delta_{ij}.
```
Normally distributed Gaussian variates are sampled using `randn()`.
## Husimi-Q representation
The [Husimi-Q representation](https://en.wikipedia.org/wiki/Husimi_Q_representation)
## Wigner representation
