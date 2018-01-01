"""
```julia
a,ā = fock(n,N;dist)
```
samples phase space distribution for a Fock state
`n` is number of fock state

`N` is number of samples

`dist` is distribution. Can be either `W` or `+P`(default)
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
    p = 0.5*sqrt(2*n+1+2*sqrt(n^2+n))
    q = 1/(4*p)
    α = (p + q*randn(N)).*exp.(2π*im*rand(N))
    ᾱ = conj(α)
    return α,ᾱ
else error("distribution not implemented")
end
end
