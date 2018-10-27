"""
    a,ā = squeezed(β,ϵ,N;dist="+P")

Sample the phase-space distribution for a squeezed state.

`β`: coherent (complex) amplitude.

`ϵ`: complex valued squeezing parameter.

`N`: number of samples.

`dist`: phase-space distribution; can be `W` or `+P`.

For standard `P,Q,W` distributions, `a` and `ā` are complex conjugate, while for `+P` etc,
`a` and `ā` are independent variables.
"""
function squeezed(β,ϵ,N;dist="+P")
if dist=="+P"
    r = abs(ϵ)
    ϕ = angle(ϵ)/2
    γ = (randn(N) .+ im*randn(N))/sqrt(2)
    ν = sqrt(exp(-r)*cosh(r)/2)*randn(N) .+ im*sqrt(exp(r)*cosh(r)/2)*randn(N)
    α = β .+ exp(im*ϕ)*ν .+ γ
    ᾱ = conj(β) .+ exp(-im*ϕ)*conj(ν) .- conj(γ)
    return α, ᾱ
elseif dist=="W"
    r = abs(ϵ)
    ϕ = angle(ϵ)/2
    α = β + 0.5*(randn(N)*exp(-r) .+ im*randn(N)*exp(r))*exp(-im*ϕ)
    ᾱ = conj(α)
    return α, ᾱ
else error("distribution not implemented")
end
end
