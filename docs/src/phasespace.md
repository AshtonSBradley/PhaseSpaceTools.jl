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
α,α⁺ = glauberP(state,N)
```
In any other distribution the coherent state is more interesting to sample.

Note here that $\alpha$ and $\alpha^+$ are simple complex conjugates. In any of the doubled-phase space representations such as the [positive-P representation](@ref), they are independent variables.


## Husimi-Q representation
The [Husimi-Q representation](https://en.wikipedia.org/wiki/Husimi_Q_representation) has the advantage that it always exists as a non-singular distribution, and is always non-negative. It is defined as the diagonal matrix element of the density matrix in the coherent state basis

```math
Q(\alpha,\alpha^*)=\frac{1}{\pi}\langle \alpha|\rho|\alpha\rangle
```

Many of its properties are very direct to prove. For example, moments of $Q$ generate anti-normally ordered operator averages:
```math
\int d^2\alpha\;(\alpha^*)^p\alpha^qQ(\alpha,\alpha^*)=\langle a^q (a^\dag)^q\rangle
```
as may be shown using properties of the trace.

In this code we make frequent use of the $Q$ function to sample other distributions. In particular, we use the fact that $P,W,Q$ generate operator averages that are normally ordered, symmetrically ordered, and anti-normally ordered respectively. This may also be expressed in the form of a convolution relationship between distributions, where $W$ and $Q$ are convolutionally broadened relative to the narrower $P$ distribution.

In the next section we outline the main application of this relationship, which is sampling particular +P distributions for states where a $Q$ function may be easily sampled. 

## Positive-P representation
In the Positive-P representation the simples way to sample the state is again as a point on the complex plane

```julia
α,α⁺ = positiveP(state,N)
```

Here ``α,α⁺`` are no longer in general complex conjugates. Even if they are initial conjugates, we emphasize that this is no longer the case in subsequent dynamical evolution of the ensemble of phase space points. Furthermore, unlike the Glauber-P distribution, the +P distribution for a particular quantum state is not unique.

In this package we primarily focus on single mode problems, taking the general approach of sampling the Husimi-Q funcion, and then exploiting a simple convolutional relationship between the Husimi-Q and a particular choice of +P representation.

A particular choice of +P that can be made for any single mode density matrix is given by
```math
P(\alpha,\alpha^+) = \frac{1}{4\pi^2}\langle (\alpha+(\alpha^+)^*)/2|\rho|(\alpha+(\alpha^+)^*)/2\rangle e^{-|\alpha-(\alpha^+)^*|^2/4}
```
We can then use the transformation
```math
\mu = (\alpha+(\alpha^+)^*)/2,\quad\quad\gamma = (\alpha-(\alpha^+)^*)/2,
```
with inverse
```math
\alpha = \mu + \gamma,\quad\quad \alpha^+ = \mu^* - \gamma^*.
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

We should emphasize that all states with $Q$ function sampling in this package are also available in the +P representation via this approach. To see details of sampling for a given distribution, see `src/methods.jl`.
## Wigner-representation
The Wigner function generates symmetrically ordered operator averages from its moments. It is often used due to being well suited for identifying semi-classical approximations. Many of the simpler Wigner functions are well known and straightforward to sample. However, when the $W$ function has negative regions, or oscillates rapidly, sampling can be more challenging. In this package we implement previously developed methods for sampling Fock states, and also develop a new method to sample Fock states in the positive Wigner representation, a doubled-phase space representation that allows the distribution to remain positive everywhere. 

## Positive Wigner representation
Just as the phase space of the P-representation can be doubled to give the positive-P, we can also double the phase space of the W-representation, to give a positive-W representation. The relationship between positive-P and positive-W is a convolution
```math
W_+(\alpha,\alpha^+)=\frac{2}{\pi}\int d^2\lambda e^{-2|\lambda|^2}P_+(\alpha+\lambda,\alpha^++\lambda^*).
```
We can make a change of variables to $\mu, \gamma$, by defining $\lambda = \mu'-(\alpha+\alpha^+)/2$, and
```math
\alpha = \mu + \gamma,\quad\quad \alpha^+ = \mu^* - \gamma^*
```
with inverse
```math
\mu = \frac{1}{2}\left(\alpha+(\alpha^+)^*\right),\quad\quad \gamma = \frac{1}{2}\left(\alpha-(\alpha^+)^*\right).
```
Using the specific form of $P_+$ in terms of the $Q$ function, we have
```math
W_+(\alpha,\alpha^+)=\frac{2}{\pi}\int d^2\mu' e^{-2|\mu'-(\alpha+(\alpha^+)^*)/2|^2}\frac{e^{-|\alpha-(\alpha^+)^*|^2/4}}{4\pi^2}Q(\mu',\mu'^*).
```
We note that
```math
\int d^2\alpha\int d^2\alpha^+ W_+(\alpha,\alpha^+)=\int d^2\gamma\int d^2\mu W_+(\alpha,\alpha^+)4
```
where the factor 4 is the Jacobian of the transformation. Hence, we define
```math
\tilde W_+(\gamma,\mu)\equiv 4W_+(\alpha(\gamma,\mu),\alpha^+(\gamma,\mu)),
```
and we thus have a particular form for the positive-W of any single mode bosonic quantum state
```math
\tilde W_+(\gamma,\mu) = \frac{e^{-|\gamma|^2}}{\pi}\frac{2}{\pi^2}\int d^2\mu' e^{-2|\mu'-\mu|^2}Q(\mu',\mu'^*)
```
together with the inverse given above
```math
W_+(\alpha,\alpha^+)=\frac{1}{4}\tilde W_+(\gamma(\alpha,\alpha^+),\mu(\alpha,\alpha^+))
```
### Fock state in positive-W
Using the Q-representation for the Fock state $\rho = |n\rangle \langle n|$, and polar coordinates $\mu=\rho e^{i\phi}$, $\mu'=r e^{i\theta}$, we can carry out the convolution
```math
\frac{2}{\pi^2}\int d^2\mu'e^{-2|\mu'|^2-2|\mu|^2+2(\mu\mu'^*+\mu^*\mu')}\frac{|\mu'|^{2n}e^{-|\mu'|^2}}{n!}
=
\frac{2}{\pi^2}\frac{e^{-2|\mu|^2}}{n!}\int_0^\infty dr r^{2n+1}e^{-3r^2}\underbrace{\int_0^{2\pi}d\theta e^{4r\rho\cos{(\theta-\phi)}}}_{2\pi I_0(4r\rho)}.
```
Using the relation between Laguerre polynomials and Bessel function
```math
\int_0^\infty dr r^{2n+1}e^{-3r^2}I_0(4r\rho)=\frac{1}{2}\frac{n!}{3^{n+1}}L_{-(n+1)}(4\rho^2/3)
```
and the identity $L_{-n}(x)=e^xL_{n-1}(-x)$, we arrive at
```math
\tilde W_+(\gamma,\mu) = \frac{e^{-|\gamma|^2}}{\pi}\frac{2e^{-2|\mu|^2/3}}{\pi 3^{n+1}}L_n(-4|\mu|^2/3)
```
The $\gamma$ distribution is easily sampled. For the distribution of $\mu\in \mathbb{C}$, as it only depends on $|\mu|\in [0,\infty)$, we can sample it via the correctly normalized distribution for $|\mu|=x$
```math
P_{|\mu|}(x)=\frac{4 x e^{-2x^2/3}}{3^{n+1}}L_n(-4x^2/3),
```
such that $\int_0^\infty \;dx P_{|\mu|}(x)=1$, 
and then apply a uniformly sampled random phase $\theta \in [0,2\pi)$ to construct samples on the complex plane
```math
\mu = x e^{i\theta}.
```
