"""
    a,ā = fock(n,N;dist=:posP)

Sample the phase-space distribution for a Fock state.

`n`: number of the fock state |n⟩.

`N`: number of samples.

`dist`: phase-space distribution can be either `:W`, `:posW`, or `:posP`.

For standard `P,Q,W` distributions, `a` and `ā` are complex conjugate, while for `+P` etc,
`a` and `ā` are independent variables.

## Wigner representation
* The standard `W` is sampled using an approximation that reproduces operator averages accuratley for large `n`, but neglects some quantum correlations.

* The `+W` sampling is carried out in a doubled phase space, where the distribution is positive semi-definite, and hence may be sampled exactly for any `n`.

* Fast evaluation is achieved for large `n (>320)` using an accurate asymptotic expansion .

# External links

___Numerical representation of quantum states in the positive-P and Wigner representations___, M. K. Olsen, A. S. Bradley, [Optics Communications 282, 3924 (2009)](http://dx.doi.org/10.1016/j.optcom.2009.06.033)

"""
function fock(n,N;dist=:posP)
if dist==:posP
    γ = randnc(N)
    d = Gamma(n+1,1)
    z = rand(d,N)
    μ = sqrt.(z).*exp.(2π*im*rand(N))
    α = μ .+ γ
    ᾱ = conj(μ .- γ)
    return α, ᾱ
elseif dist==:W
    n<10 && warn("Fock state sampling for W is only valid for n ≫ 1.")
    p = 0.5*sqrt(2*n+1+2*sqrt(n^2+n))
    q = 1/(4*p)
    α = (p .+ q*randn(N)).*exp.(2π*im*rand(N))
    ᾱ = conj(α)
    return α, ᾱ
elseif dist==:posW
    if n<=320
    γ = randnc(N)
    x1 = max(0,sqrt(n)-5); x2 = sqrt(n)+5
    (n==0||n==1) ? Pmax=0.71 : Pmax=0.6
    z = reject(x->plaguerre.(x,n),[x1,x2],N,Pmax)
    μ = z.*exp.(2π*im*rand(N))
    α = μ .+ γ
    ᾱ = conj(μ .- γ)
    return α, ᾱ
    else
    γ = randnc(N)
    x1= max(0,sqrt(n)-5); x2 = sqrt(n)+5
    z = reject(x->plaguerre_asymptotic.(x,n),[x1,x2],N,0.6)
    μ = z.*exp.(2π*im*rand(N))
    α = μ .+ γ
    ᾱ = conj(μ .- γ)
    return α, ᾱ
    end
else error("distribution not implemented")
end
end

"""
    plaguerre(x,n)

Define the radial phase-space probability distribution for sampling a fock state in the
`+W` representation, where the probability density is
```math
P(x,n) \\equiv \\frac{4}{3^{n+1}}\\exp{(-2x^2/3)}xL_n(-4x^2/3)
```

`x`: location to evaluate probability.

`n`: the number of the fock state.

"""
plaguerre(x,n) = (4/3)*exp(-2*x^2/3 + log(x) + log(laguerren(-4*x^2/3,n)) - n*log(3))

"""
    plaguerre_asymptotic(x,n)

Define the radial phase-space probability distribution for sampling a fock state in the
`+W` representation. In the large `n` limit the distribution is sampled using the ansymptotic expansion
```math
P(x,n)=\\frac{1}{\\sqrt{\\pi}}\\exp{(-[x-\\sqrt{n+1}]^2)}
```
where
`x`: point to evaluate probability.

`n`: number of the fock state |n⟩.

"""
plaguerre_asymptotic(x,n) = exp(-(x-sqrt(n+1))^2)/sqrt(pi)

#Evaluate Laguerre polynomial by calling GSL
laguerren(x,n)=sf_laguerre_n(n,0.0,x)
