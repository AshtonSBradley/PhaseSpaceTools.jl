"""
    a,ā = wigner(state <: State,N)

Generate `N` samples from wigner phase-space distribution for `state`.

Moments of the Wigner distribution generate symmetrically ordered quantum operator averages.
"""
function wigner(state::Coherent,N)
    @unpack β = state
    α = β .+ randnc(N)/sqrt(2)
    ᾱ = conj(α)
    return α, ᾱ
end

"""
    a,ā = husimiQ(state <: State,N)

Generate `N` samples from the Husimi-Q phase-space distribution for `state`.

Moments of the Husimi-Q distribution generate quantum operator averages that are anti-normally ordered.
"""
function husimiQ(state::Coherent,N)
    @unpack β = state
    α = β .+ randnc(N)
    ᾱ = conj(α)
    return α, ᾱ
end

"""
    a,ā = positiveP(state <: State,N)

Generate `N` samples from the positive-P phase-space distribution for `state`.

Moments of the positive-P distribution generate quantum operator averages that are normally ordered.
"""
function positiveP(state::T,N) where T
    μ,μ̄ = husimiQ(state,N)
    γ = randnc(N)
    α = μ .+ γ
    ᾱ = conj(μ .- γ)
    return α, ᾱ
end

function positiveP(state::Coherent,N)
    @unpack β = state
    α = β*ones(N)
    ᾱ = conj(α)
    return α,ᾱ
end

"""
    a,ā = glauberP(state <: State,N)

Generate `N` samples from the Glauber-Sudarshan-P phase-space distribution for `state`.

Moments of the Glauber-Sudarshan-P distribution generate quantum operator averages that are normally ordered.
"""
glauberP(state::Coherent,N) = positiveP(state,N)

function glauberP(state::Thermal,N)
    @unpack β,n̄ = state
    α = β .+ sqrt(n̄)*randnc(N)
    ᾱ = conj(α)
    return α, ᾱ
end

function wigner(state::Thermal,N)
    @unpack β,n̄ = state
    α = β .+ sqrt(n̄+.5)*randnc(N)
    ᾱ = conj(α)
    return α, ᾱ
end

function husimiQ(state::Thermal,N)
    @unpack β,n̄ = state
    α = β .+ sqrt(n̄+1.0)*randnc(N)
    ᾱ = conj(α)
    return α, ᾱ
end

function husimiQ(state::Squeezed,N)
    @unpack β,ϵ = state
    r = abs(ϵ)
    ϕ = angle(ϵ)/2
    ν = sqrt(exp(-r)*cosh(r)/2)*randn(N) .+ im*sqrt(exp(r)*cosh(r)/2)*randn(N)
    α = β .+ exp(im*ϕ)*ν
    ᾱ = conj(α)
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

function husimiQ(state::Crescent,N)
    @unpack β,ϵ,q = state
    r = abs(ϵ)
    ϕ = angle(ϵ)/2
    α = β .+ (randn(N)*exp(-r) .+ im*randn(N)*exp(r))*exp(-im*ϕ)/sqrt(2)
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

function husimiQ(state::Fock,N)
    @unpack n = state
    d = Gamma(n+1,1)
    z = rand(d,N)
    α = sqrt.(z).*exp.(2π*im*rand(N))
    ᾱ = conj(α)
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

"""
    a,ā = positiveW(state <: State,N)

Generate `N` samples from the positive-W phase-space distribution for `state`.

Implemented states are

- `Fock(N)`

The moments of the positive-W distribution generate quantum operator averages that are symmmetrically ordered.
"""
function positiveW(state::Fock,N)
    @unpack n = state
    if n<=320
    x1 = max(0,sqrt(n)-5); x2 = sqrt(n)+5
    (n==0||n==1) ? Pmax=0.71 : Pmax=0.6
    z = reject(x->plaguerre.(x,n),[x1,x2],N,Pmax)
    μ = z.*exp.(2π*im*rand(N))
    γ = randnc(N)
    α = μ .+ γ
    ᾱ = conj(μ .- γ)
    return α, ᾱ
    else
    x1 = max(0,sqrt(n)-5); x2 = sqrt(n)+5
    z = reject(x->plaguerre_asymptotic.(x,n),[x1,x2],N,0.6)
    μ = z.*exp.(2π*im*rand(N))
    γ = randnc(N)
    α = μ .+ γ
    ᾱ = conj(μ .- γ)
    return α, ᾱ
end
end

function glauberP(state::Bogoliubov,N)
    @unpack u,v,n̄ = state
    b,b̄ = glauberP(Thermal(0.0,n̄),N)
    a = u*b + conj(v)*b̄
    ā = conj.(a)
    return a,ā
end

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

# #TODO
# function positiveW(state::Coherent,N)
#     @unpack β = state
#     α = β .+ randnc(N)/sqrt(2)
#     ᾱ = conj(α)
#     return α, ᾱ
# end
