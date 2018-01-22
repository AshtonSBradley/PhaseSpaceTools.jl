"""
    a,ā = fock(n,N;dist="+P")

Sample the phase-space distribution for a Fock state.

`n`: number of the fock state |n⟩.

`N`: number of samples.

`dist`: phase-space distribution can be either `W`, `+W`, or `+P`.

For standard `P,Q,W` distributions, `a` and `ā` are complex conjugate, while for `+P` etc,
`a` and `ā` are independent variables.

## Wigner representation
* The standard `W` is sampled using an approximation that reproduces operator averages accuratley for large `n`, but neglects some quantum correlations.

* The `+W` sampling is carried out in a doubled phase space, where the distribution is positive semi-definite, and hence may be sampled exactly for any `n`. Fast evaluation is achieved for large `n` using an accurate asymptotic expansion `(n>320)`.

External links

[Numerical representation of quantum states in the positive-P and Wigner representations, Olsen, Bradley, Optics Communications 282, 3924 (2009)](http://dx.doi.org/10.1016/j.optcom.2009.06.033)

"""
function fock(n,N;dist="+P")
if dist=="+P"
    γ = (randn(N)+im*randn(N))/sqrt(2)
    d = Gamma(n+1,1)
    z = rand(d,N)
    μ = sqrt.(z).*exp.(2π*im*rand(N))
    α = μ + γ
    ᾱ = conj(μ - γ)
    return α, ᾱ
elseif dist=="W"
    n<10 && warn("Fock state sampling for W is only valid for n ≫ 1.")
    p = 0.5*sqrt(2*n+1+2*sqrt(n^2+n))
    q = 1/(4*p)
    α = (p + q*randn(N)).*exp.(2π*im*rand(N))
    ᾱ = conj(α)
    return α, ᾱ
elseif dist=="+W"
    if n<=320
    γ = (randn(N)+im*randn(N))/sqrt(2)
    x1= max(0,sqrt(n)-5); x2 = sqrt(n)+5
    (n==0||n==1) ? Pmax=0.71 : Pmax=0.6
    z = reject(x->plaguerre.(x,n),[x1,x2],N,Pmax)
    μ = z.*exp.(2π*im*rand(N))
    α = μ + γ
    ᾱ = conj(μ - γ)
    return α, ᾱ
    else
    γ = (randn(N)+im*randn(N))/sqrt(2)
    x1= max(0,sqrt(n)-5); x2 = sqrt(n)+5
    z = reject(x->plaguerre_asymptotic.(x,n),[x1,x2],N,0.6)
    μ = z.*exp.(2π*im*rand(N))
    α = μ + γ
    ᾱ = conj(μ - γ)
    return α, ᾱ
    end
else error("distribution not implemented")
end
end
