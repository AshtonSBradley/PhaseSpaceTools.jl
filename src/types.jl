"""
    State

Abstract supertype for all sampled states: `state <: State`.

# Examples
## Find all states that may be sampled

```julia-repl
julia> subtypes(State)

7-element Vector{Any}:
Bogoliubov
Coherent
Crescent
Fock
Squeezed
SqueezedTwoMode
Thermal
```

## Create and sample a particular state (vacuum)

```julia-repl 
julia> s = Fock(0)
Fock(0)
julia> wigner(s,100)
┌ Warning: Fock state sampling for W is only valid for n ≫ 1.
```

Here a warning is generated since Fock sampling is not well defined for small `n`. 
    
A simpler way to sample the vacuum is 

```julia-repl 
julia> s = Coherent(0) 
Coherent(0.0 + 0.0im)  # type conversion to ComplexF64.

julia> wigner(s,100)
(ComplexF64[0.33820868828162637 + 0.4407579103538181im, 0.057183146091823775 - 0.2772571883006981im, ...
```
generating two vectors of sampled points `α,α⁺` in the complex plane. In this case, `α = conj(α⁺)`, as we are not working with a doubled phase space.
"""    
abstract type State end

"""
    Coherent(β)

Create a coherent state with complex amplitude `β ∈ ℂ`.
"""
struct Coherent <: State
    β::Complex{Float64}
end

Coherent(β::Number) = Coherent(ComplexF64(β))

"""
    Fock(n)

Create a Fock state for particle number `n ∈ ℕ₀`.
"""
struct Fock <: State
    n::Int64
    function Fock(n::Int64)
        n >= 0 || throw(ArgumentError("Fock state requires n >= 0."))
        return new(n)
    end
end

Fock(n::Integer) = Fock(Int64(n))

"""
    Crescent(β,ϵ,q)

Create a Crescent state with parameters `β`, `ϵ`, `q`.
"""
struct Crescent <: State
    β::Complex{Float64}
    ϵ::Complex{Float64}
    q::Float64
end

Crescent(β::Number, ϵ::Number, q::Real) = Crescent(ComplexF64(β), ComplexF64(ϵ), Float64(q))

"""
    Squeezed(β,ϵ)

Create a Squeezed state with parameters `β`, `ϵ`.
"""
struct Squeezed <: State
    β::Complex{Float64}
    ϵ::Complex{Float64}
end

Squeezed(β::Number, ϵ::Number) = Squeezed(ComplexF64(β), ComplexF64(ϵ))

"""
    SqueezedTwoMode(r,ϕ)

Create a two-mode squeezed state with parameters `r`, `ϕ`.
"""
struct SqueezedTwoMode <: State
    r::Float64
    ϕ::Float64
    function SqueezedTwoMode(r::Float64, ϕ::Float64)
        r >= 0 || throw(ArgumentError("SqueezedTwoMode requires r >= 0."))
        return new(r, ϕ)
    end
end

SqueezedTwoMode(r::Real, ϕ::Real) = SqueezedTwoMode(Float64(r), Float64(ϕ))

"""
    Thermal(β,n̄)

Create a Thermal state with parameters `β`, `n̄`.
"""
struct Thermal <: State
    β::Complex{Float64}
    n̄::Float64
    function Thermal(β::Complex{Float64}, n̄::Float64)
        n̄ >= 0 || throw(ArgumentError("Thermal state requires n̄ >= 0."))
        return new(β, n̄)
    end
end

Thermal(β::Number, n̄::Real) = Thermal(ComplexF64(β), Float64(n̄))

"""
    Bogoliubov(u,v,n̄)

Create a Bogoliubov state with parameters `u`, `v`, `n̄`.
"""
struct Bogoliubov <: State
    u::Complex{Float64} 
    v::Complex{Float64}
    n̄::Float64
    function Bogoliubov(u::Complex{Float64}, v::Complex{Float64}, n̄::Float64)
        n̄ >= 0 || throw(ArgumentError("Bogoliubov state requires n̄ >= 0."))
        isapprox(abs2(u) - abs2(v), 1.0; atol=1e-8, rtol=1e-8) ||
            throw(ArgumentError("Bogoliubov state requires |u|^2 - |v|^2 ≈ 1."))
        return new(u, v, n̄)
    end
end

Bogoliubov(u::Number, v::Number, n̄::Real) = Bogoliubov(ComplexF64(u), ComplexF64(v), Float64(n̄))
