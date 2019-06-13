# old methods

function positiveP(state::Thermal,N)
    μ,μ̄ = husimiQ(state,N)
    γ = crandn(N)
    α = μ .+ γ
    ᾱ = conj(μ .- γ)
    return α, ᾱ
end

function positiveP(state::Squeezed,N)
    μ,μ̄ = husimiQ(state,N)
    γ = crandn(N)
    α = μ .+ γ
    ᾱ = conj(μ .- γ)
    return α, ᾱ
end

function positiveP(state::Crescent,N)
    μ,μ̄ = husimiQ(state,N)
    γ = crandn(N)
    α = μ .+ γ
    ᾱ = conj(μ .- γ)
    return α,ᾱ
end

function positiveP(state::Fock,N)
    μ,μ̄ = husimiQ(state,N)
    γ = crandn(N)
    α = μ .+ γ
    ᾱ = conj(μ .- γ)
    return α, ᾱ
end

function positiveP(state::Bogoliubov,N)
    μ,μ̄ = husimiQ(state,N)
    γ = crandn(N)
    α = μ .+ γ
    ᾱ = conj(μ .- γ)
    return α, ᾱ
end
