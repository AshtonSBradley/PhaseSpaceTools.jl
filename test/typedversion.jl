using Distributions, Statistics
import GSL:sf_laguerre_n

abstract type State end

struct Coherent <: State
    β::Complex{Float64}
end

struct Fock <: State
    n::Int
end

struct Crescent <: State
    β::Complex{Float64}
    ϵ::Float64
    q::Float64
end

struct Squeezed <: State
    β::Complex{Float64}
    ϵ::Float64
end

struct Thermal <: State
    β::Complex{Float64}
    n̄::Float64
end

# helpers
function randnc(args...)
    return randn(ComplexF64,args...)
end

function reject(P,w,N,Pmax)
    a,b = w
    samples = Array{Float64}(undef,0)
    sizehint!(samples,N)
    while length(samples) < N
        y = a + rand()*(b - a)
        z = rand()*Pmax
        z < P(y) && push!(samples,y)
    end
    return samples
end

plaguerre(x,n) = (4/3)*exp(-2*x^2/3 + log(x) + log(laguerren(-4*x^2/3,n)) - n*log(3))
plaguerre_asymptotic(x,n) = exp(-(x-sqrt(n+1))^2)/sqrt(pi)
laguerren(x,n) = sf_laguerre_n(n,0.0,x)

#sampling methods
function positiveP(state::Thermal,N)
    β,n̄ = state.β,state.n̄
    α = β .+ sqrt(n̄)*randnc(N)
    ᾱ = conj(α)
    return α, ᾱ
end

glauberP(state::Thermal,N) = positiveP(state,N)

function husimiQ(state::Thermal,N)
    α = β .+ sqrt(n̄+1.0)*randnc(N)
    ᾱ = conj(α)
    return α, ᾱ
end

function wigner(state::Thermal,N)
    β,n̄ = state.β,state.n̄
    α = β .+ sqrt(n̄+.5)*randnc(N)
    ᾱ = conj(α)
    return α, ᾱ
end

function positiveP(state::Squeezed,N)
    β,ϵ = state.β,state.ϵ
    r = abs(ϵ)
    ϕ = angle(ϵ)/2
    γ = randnc(N)
    ν = sqrt(exp(-r)*cosh(r)/2)*randn(N) .+ im*sqrt(exp(r)*cosh(r)/2)*randn(N)
    α = β .+ exp(im*ϕ)*ν .+ γ
    ᾱ = conj(β) .+ exp(-im*ϕ)*conj(ν) .- conj(γ)
    return α, ᾱ
end

function wigner(state::Squeezed,N)
    r = abs(ϵ)
    ϕ = angle(ϵ)/2
    α = β .+ 0.5*(randn(N)*exp(-r) .+ im*randn(N)*exp(r))*exp(-im*ϕ)
    ᾱ = conj(α)
    return α, ᾱ
end

function positiveP(state::Coherent,N)
    β = state.β
    α = β*ones(N)
    ᾱ = conj(α)
    return α,ᾱ
end

glauberP(state::Coherent,N) = positiveP(state,N)

function positiveW(state::Coherent,N)
    β = state.β
    α = β .+ randnc(N)/sqrt(2)
    ᾱ = conj(α)
    return α, ᾱ
end

wigner(state::Coherent,N) = positiveW(state::Coherent,N)

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

function wigner(state::Fock,N)
    n = state.n
    n < 10 && warn("Fock state sampling for W is only valid for n ≫ 1.")
    p = 0.5*sqrt(2*n+1+2*sqrt(n^2+n))
    q = 1/(4*p)
    α = (p .+ q*randn(N)).*exp.(2π*im*rand(N))
    ᾱ = conj(α)
    return α, ᾱ
end

function positiveW(state::Fock,N)
    n = state.n
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
end

# tests

N = 1000000
β = 100*randnc()
state = Coherent(β)
a,ā = wigner(state,N)
meanb = mean(a)
n̄ = real(mean(@. ā*a)-.5)
Vn = mean(@. ā^2*a^2)-mean(@. a*ā)-n̄^2 |> real

state = Fock(15)
α,ᾱ = positiveP(state,100)
@time α,ᾱ = wigner(state,1000)
@time α,ᾱ = positiveW(state,1000000)

state = Squeezed(20. +im*1.9,2.0)



methods(positiveP)
