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
    α = β .+ crandn(N)/sqrt(2)
    ᾱ = conj(α)
    return α, ᾱ
end

wigner(state::Coherent,N) = positiveW(state::Coherent,N)

function positiveP(state::Thermal,N)
    @unpack β,n̄ = state
    α = β .+ sqrt(n̄)*crandn(N)
    ᾱ = conj(α)
    return α, ᾱ
end

glauberP(state::Thermal,N) = positiveP(state,N)

function husimiQ(state::Thermal,N)
    @unpack β,n̄ = state
    α = β .+ sqrt(n̄+1.0)*crandn(N)
    ᾱ = conj(α)
    return α, ᾱ
end

function wigner(state::Thermal,N)
    @unpack β,n̄ = state
    α = β .+ sqrt(n̄+.5)*crandn(N)
    ᾱ = conj(α)
    return α, ᾱ
end

function positiveP(state::Squeezed,N)
    @unpack β,ϵ = state
    r = abs(ϵ)
    ϕ = angle(ϵ)/2
    γ = crandn(N)
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
    γ = crandn(N)
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
    γ = crandn(N)
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
    x1 = max(0,sqrt(n)-5); x2 = sqrt(n)+5
    z = reject(x->plaguerre_asymptotic.(x,n),[x1,x2],N,0.6)
    μ = z.*exp.(2π*im*rand(N))
    α = μ .+ γ
    ᾱ = conj(μ .- γ)
    return α, ᾱ
end
end

function positiveP(state::Bogoliubov,N)
    @unpack u,v,n̄ = state
    b,b̄ = positiveP(Thermal(0.0,n̄),N)
    a = u*b + conj(v)*b̄
    ā = conj.(a)
    return a,ā
end

glauberP(state::Bogoliubov,N) = positiveP(state,N)

function wigner(state::Bogoliubov,N)
    @unpack u,v,n̄ = state
    b,b̄ = wigner(Thermal(0.0,n̄),N)
    a = u*b + conj(v)*b̄
    ā = conj.(a)
    return a,ā
end

function husimiQ(state::Bogoliubov,N)
    @unpack u,v,n̄ = state
    b,b̄ = husimiQ(Thermal(0.0,n̄),N)
    a = u*b + conj(v)*b̄
    ā = conj.(a)
    return a,ā
end
