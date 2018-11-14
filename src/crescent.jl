"""
    a,ā = crescent(β,ϵ,q,N;dist=:posP)

Sample the phase-space distribution for a crescent state.
First samples a squeezed state, then introduces a searing factor
in phase-space in the form of a gaussian distributed random phase.

`β` is coherent (complex) amplitude.

`ϵ`: squeezing paramter.

`q`: shearing parameter.

`N`: number of samples.

`dist`: phase-space distribution. Can be `:posP`,`:Q` or `:W`.

For standard `P,Q,W` distributions, `a` and `ā` are complex conjugate, while for `+P` etc,
`a` and `ā` are independent variables.
"""
function crescent(β,ϵ,q,N;dist=:posP)
if dist==:posP
    r = abs(ϵ)
    ϕ = angle(ϵ)/2
    μ = β .+ (randn(N)*exp(-r)+im*randn(N)*exp(r))/sqrt(2)
    μ .*= exp.(im*q*randn(N))
    γ = crandn(N)
    α = μ .+ γ
    ᾱ = conj(μ .- γ)
    return α,ᾱ
elseif dist==:Q
    r = abs(ϵ)
    ϕ = angle(ϵ)/2
    α = β .+ (randn(N)*exp(-r) .+ im*randn(N)*exp(r))/sqrt(2)
    α = α.*exp.(im*q*randn(N))
    ᾱ = conj(α)
    return α,ᾱ
elseif dist==:W
    r = abs(ϵ)
    ϕ = angle(ϵ)/2
    α = β .+ 0.5*(randn(N)*exp(-r) .+ im*randn(N)*exp(r))*exp(-im*ϕ)
    α = α.*exp.(im*q*randn(N))
    ᾱ = conj(α)
    return α, ᾱ
else error("distribution not implemented")
end
end
