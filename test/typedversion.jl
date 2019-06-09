using Distributions, Statistics

abstract type State end

struct Coherent <: State
    β::Complex{Float64}
end

struct Fock <: State
    n::Int
end

function randnc(args...)
    return randn(ComplexF64,args...)
end

function positiveP(state::Coherent,N)
    β = state.β
    α = β*ones(N)
    ᾱ = conj(α)
    return α,ᾱ
end

function wigner(state::Coherent,N)
    β = state.β
    α = β .+ randnc(N)/sqrt(2)
    ᾱ = conj(α)
    return α,ᾱ
end

function positiveW(state::Coherent,N)
    β = state.β
    α = β .+ randnc(N)/sqrt(2)
    ᾱ = conj(α)
    return α, ᾱ
end

# test coherent

state = Coherent(10.2+im)

α,ᾱ = wigner(state,1000)

struct Crescent <: State
    β::Complex{Float64}
    ϵ::Float64
    q::Float64
end

function positiveP(state::Crescent,N)
    β,ϵ,q = state.β,state.ϵ,state.q
    r = abs(ϵ)
    ϕ = angle(ϵ)/2
    μ = β .+ (randn(N)*exp(-r)+im*randn(N)*exp(r))/sqrt(2)
    μ .*= exp.(im*q*randn(N))
    γ = randnc(N)
    α = μ .+ γ
    ᾱ = conj(μ .- γ)
    return α,ᾱ
end

function husimiQ(state::Crescent,N)
    β,ϵ,q = state.β,state.ϵ,state.q
    r = abs(ϵ)
    ϕ = angle(ϵ)/2
    α = β .+ (randn(N)*exp(-r) .+ im*randn(N)*exp(r))/sqrt(2)
    α = α.*exp.(im*q*randn(N))
    ᾱ = conj(α)
    return α,ᾱ
end

function wigner(state::Crescent,N)
    β,ϵ,q = state.β,state.ϵ,state.q
    r = abs(ϵ)
    ϕ = angle(ϵ)/2
    α = β .+ 0.5*(randn(N)*exp(-r) .+ im*randn(N)*exp(r))*exp(-im*ϕ)
    α = α.*exp.(im*q*randn(N))
    ᾱ = conj(α)
    return α, ᾱ
end

# test fock

state = Fock(12)

function positiveP(state::Fock,N)
    n = state.n
    γ = randnc(N)
    d = Gamma(n+1,1)
    z = rand(d,N)
    μ = sqrt.(z).*exp.(2π*im*rand(N))
    α = μ .+ γ
    ᾱ = conj(μ .- γ)
    return α, ᾱ
end

α,ᾱ = positiveP(state,100)

function wigner(state::Fock,N)
    n = state.n
    n < 10 && warn("Fock state sampling for W is only valid for n ≫ 1.")
    p = 0.5*sqrt(2*n+1+2*sqrt(n^2+n))
    q = 1/(4*p)
    α = (p .+ q*randn(N)).*exp.(2π*im*rand(N))
    ᾱ = conj(α)
    return α, ᾱ
end

α,ᾱ = wigner(state,100)

function fock(n,N;dist=:posP)
if dist==:posP
    γ = crandn(N)
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
    γ = crandn(N)
    x1 = max(0,sqrt(n)-5); x2 = sqrt(n)+5
    (n==0||n==1) ? Pmax=0.71 : Pmax=0.6
    z = reject(x->plaguerre.(x,n),[x1,x2],N,Pmax)
    μ = z.*exp.(2π*im*rand(N))
    α = μ .+ γ
    ᾱ = conj(μ .- γ)
    return α, ᾱ
    else
    γ = crandn(N)
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
