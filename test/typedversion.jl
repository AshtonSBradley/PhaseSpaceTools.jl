using Distributions, Statistics, Parameters
import GSL:sf_laguerre_n

# helpers
randnc(args...) = randn(ComplexF64,args...)

abstract type State end

struct Coherent <: State
    β::Complex{Float64}
end

struct Fock <: State
    n::Int64
end

struct Crescent <: State
    β::Complex{Float64}
    ϵ::Float64
    q::Float64
end

struct Squeezed <: State
    β::Complex{Float64}
    ϵ::Complex{Float64}
end

struct Thermal <: State
    β::Complex{Float64}
    n̄::Float64
end

struct Bogoliubov <: State
    u::Complex{Float64}
    v::Complex{Float64}
    n̄::Float64
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
function positiveP(state::Coherent,N)
    @unpack β = state
    α = β*ones(N)
    ᾱ = conj(α)
    return α,ᾱ
end

glauberP(state::Coherent,N) = positiveP(state,N)

function positiveW(state::Coherent,N)
    @unpack β = state
    α = β .+ randnc(N)/sqrt(2)
    ᾱ = conj(α)
    return α, ᾱ
end

wigner(state::Coherent,N) = positiveW(state::Coherent,N)

function positiveP(state::Thermal,N)
    @unpack β,n̄ = state
    α = β .+ sqrt(n̄)*randnc(N)
    ᾱ = conj(α)
    return α, ᾱ
end

glauberP(state::Thermal,N) = positiveP(state,N)

function husimiQ(state::Thermal,N)
    @unpack β,n̄ = state
    α = β .+ sqrt(n̄+1.0)*randnc(N)
    ᾱ = conj(α)
    return α, ᾱ
end

function wigner(state::Thermal,N)
    @unpack β,n̄ = state
    α = β .+ sqrt(n̄+.5)*randnc(N)
    ᾱ = conj(α)
    return α, ᾱ
end

function positiveP(state::Squeezed,N)
    @unpack β,ϵ = state
    r = abs(ϵ)
    ϕ = angle(ϵ)/2
    γ = randnc(N)
    ν = sqrt(exp(-r)*cosh(r)/2)*randn(N) .+ im*sqrt(exp(r)*cosh(r)/2)*randn(N)
    α = β .+ exp(im*ϕ)*ν .+ γ
    ᾱ = conj(β) .+ exp(-im*ϕ)*conj(ν) .- conj(γ)
    return α, ᾱ
end

function wigner(state::Squeezed,N)
    @unpack β,ϵ = state
    r = abs(ϵ)
    ϕ = angle(ϵ)/2
    α = β .+ 0.5*(randn(N)*exp(-r) .+ im*randn(N)*exp(r))*exp(-im*ϕ)
    ᾱ = conj(α)
    return α, ᾱ
end

function positiveP(state::Crescent,N)
    @unpack β,ϵ,q = state
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
    @unpack β,ϵ,q = state
    r = abs(ϵ)
    ϕ = angle(ϵ)/2
    α = β .+ (randn(N)*exp(-r) .+ im*randn(N)*exp(r))/sqrt(2)
    α = α.*exp.(im*q*randn(N))
    ᾱ = conj(α)
    return α,ᾱ
end

function wigner(state::Crescent,N)
    @unpack β,ϵ,q = state
    r = abs(ϵ)
    ϕ = angle(ϵ)/2
    α = β .+ 0.5*(randn(N)*exp(-r) .+ im*randn(N)*exp(r))*exp(-im*ϕ)
    α = α.*exp.(im*q*randn(N))
    ᾱ = conj(α)
    return α, ᾱ
end

function positiveP(state::Fock,N)
    @unpack n = state
    γ = randnc(N)
    d = Gamma(n+1,1)
    z = rand(d,N)
    μ = sqrt.(z).*exp.(2π*im*rand(N))
    α = μ .+ γ
    ᾱ = conj(μ .- γ)
    return α, ᾱ
end

function wigner(state::Fock,N)
    @unpack n = state
    n < 10 && warn("Fock state sampling for W is only valid for n ≫ 1.")
    p = 0.5*sqrt(2*n+1+2*sqrt(n^2+n))
    q = 1/(4*p)
    α = (p .+ q*randn(N)).*exp.(2π*im*rand(N))
    ᾱ = conj(α)
    return α, ᾱ
end

function positiveW(state::Fock,N)
    @unpack n = state
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


# test coherent
using Test, LinearAlgebra

# sample
N = 100000
α = 100*randnc()
state = Coherent(α)

a,ā = wigner(state,N)
meana = mean(a)
n̄ = real(mean(@. ā*a)-.5)
Vn = mean(@. ā^2*a^2)-mean(@. a*ā)-n̄^2 |> abs

#test
@test norm(meana - α)/abs(α) < 0.01
@test abs(n̄ - abs(α).^2)/abs(α)^2 < 0.01
@test abs(Vn - abs(α)^2)/abs(α)^2 < 0.05

#sample
n = 100
N = 100000
state = Fock(n)
a,ā = wigner(state,N)

meana = mean(a)
absa = abs(meana)
n̄ = mean(abs2.(a))-.5
Vn= abs(mean(abs2.(a).^2)-mean(abs2.(a))-n̄.^2)
rel_num_var = sqrt(abs(Vn))/abs(n̄);

#test
@test absa < 0.1
@test n̄ - n < 0.01
@test Vn < 0.01
@test rel_num_var < 0.001

#sample
n = 99
N = 100000
state = Fock(n)
a,ā = positiveW(state,N)

meana = mean(a)
absa = abs(meana)
n̄ = real(mean(a.*ā)-.5)
Vn = abs(mean(a.*a.*ā.*ā)-abs(mean(a.*ā))-n̄.^2)
rel_num_var = sqrt(abs(Vn))/abs(n̄);

#test
@test absa < 0.1
@test n̄ - n < 1
@test Vn < 50
@test rel_num_var < 0.1

#sample
n = 101
state = Fock(n)
N = 100000

a,ā = positiveW(state,N)

av_a = mean(a)
absa = abs(av_a)
n̄ = real(mean(a.*ā)-.5)
Vn= abs(mean(a.*a.*ā.*ā)-mean(a.*ā)-n̄.^2)
rel_num_var = sqrt(abs(Vn))/abs(n̄);

#test
@test absa < 0.1
@test (n̄ - n)/n < 0.01
@test Vn < 30
@test rel_num_var < 0.1

#sample
β = 10
ϕ = π/16
r = 2
ϵ = r*exp(2*im*ϕ)
state = Squeezed(β,ϵ)
N = 10000

a,ā = positiveP(state,N)

av_a = mean(a);absa = abs(av_a)
n̄ = real(mean(a.*ā))
nbar = sinh(abs(ϵ)).^2 + abs2(β)

#test
@test norm(av_a - β)/abs(β) < 0.01
@test abs(n̄ - nbar)/abs(β)^2 < 0.01

#sample
β = 10
ϕ = π/16
r = 1.5
ϵ = r*exp(2*im*ϕ)
state = Squeezed(β,ϵ)
N = 10000
a,ā = wigner(state,N)

av_a = mean(a);absa = abs(av_a)
n̄ = real(mean(a.*ā)-.5)
nbar = sinh(abs(ϵ)).^2+abs2(β)

#test
@test norm(av_a - β)/abs(β) < 0.01
@test abs(n̄ - nbar)/abs(β)^2 < 0.01



#TODO Bogoliubov
function wigner(state::Bogoliubov,N)
    @unpack u,v,n̄ = state
    b,b̄ = wigner(Thermal(0.0,n̄),N)
    a = u*b + conj(v)*b̄
    ā = conj.(a)
    return a,ā
end

#sample
N = 1000000
n̄ = 10

a = randnc()
b = randnc()
u = a+b
v = a-b
nrm = abs2(u)-abs2(v)
u /= sqrt(nrm)
v /= sqrt(nrm)

abs2(u)-abs2(v) ≈ 1.0
state = Bogoliubov(u,v,n̄)
a,ā = wigner(state,N)

# thermal mode population (zero coherent amplitude)
N̄ = real(mean(a.*ā)-.5)
(abs2(u)+abs2(v))/2

# analytic form for Bogoliubov state:

#test
@test isapprox(N̄,n̄,rtol=1e-2)

#TODO Crescent tests
