#helper functions
"""
    a = randnc(N...)

Return an array of dimension `length(N)`, containing samples of complex random
variables with mean zero and variance one:
```math
⟨|a|^2⟩ = 1; ⟨a^2⟩ = ⟨(a^*) ^2⟩ = 0,
```
Useful for creating a range of noises that show up in phase-space simulations.

Usage:

`a = randnc(10)` returns a `10-element Array{Complex{Float64},1}`.

`a = randnc(50,100)` returns a `50x100-element Array{Complex{Float64},2}`.
"""
randnc(args...) = randn(ComplexF64,args...)


"""
    x = reject(P,w,N,Pmax)

Generate `x` distributed according to a probability distribution by rejection sampling over a finite window.

`P`: probability distribution, normalized to 1.

`w=[w1,w2]`: window for sampling `x`.

`N`: number of samples.

`Pmax`: numerical upper bound for `P`: `Pmax ≧ max(P(w))`.

"""
function reject(P,w,N,Pmax)
    a,b = w
    N >= 0 || throw(ArgumentError("reject requires N >= 0."))
    a < b || throw(ArgumentError("reject requires a finite window with w[1] < w[2]."))
    isfinite(Pmax) && Pmax > 0 || throw(ArgumentError("reject requires finite Pmax > 0."))
    samples = Array{Float64}(undef,0)
    sizehint!(samples,N)
    attempts = 0
    max_attempts = max(1_000, 1_000 * max(N, 1))
    while length(samples) < N
        y = a + rand()*(b - a)
        z = rand()*Pmax
        py = P(y)
        isfinite(py) || throw(ArgumentError("reject requires P(y) to be finite over the sampling window."))
        py >= 0 || throw(ArgumentError("reject requires P(y) >= 0 over the sampling window."))
        py <= Pmax || throw(ArgumentError("reject requires Pmax >= P(y) over the sampling window."))
        z < py && push!(samples,y)
        attempts += 1
        attempts <= max_attempts || throw(ArgumentError("reject exceeded $max_attempts attempts without collecting $N samples; check the sampling window and Pmax."))
    end
    return samples
end

# Laguerre methods for +W
plaguerre(x,n) = (4/3)*exp(-2*x^2/3 + log(x) + log(laguerren(-4*x^2/3,n)) - n*log(3))
plaguerre_asymptotic(x,n) = exp(-(x-sqrt(n+1))^2)/sqrt(pi)
laguerren(x,n) = sf_laguerre_n(n,0.0,x)


"""
Complex variable SDE's in `DifferentialEquations.jl`
have complex noises by default dispatch on complex fields.
This method defines a real noise for use with the library.
"""
function realnoise(rand_vec,W,dt,rng)
for i in eachindex(rand_vec)
    rand_vec[i] = randn(rng)
    end
    rand_vec .*= sqrt(abs(dt))
end

"""
Complex variable SDE's in `DifferentialEquations.jl`
have complex noises by default dispatch on complex fields.
This method defines a real brownian bridge for use with the library.
"""
function realbridge(rand_vec,W,W0,Wh,q,h,rng)
for i in eachindex(rand_vec)
    rand_vec[i] = randn(rng)
    end
    rand_vec .= sqrt((1 .- q).*q.*abs(h)).*rand_vec.+q.*Wh
end
