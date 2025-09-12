"""
    State

Abstract supertype for all sampled states: `state <: State`.

# Examples
## Find all states that may be sampled

```julia-repl
julia> subtypes(State)

6-element Vector{Any}:
Bogoliubov
Coherent
Crescent
Fock
Squeezed
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

"""
    Fock(n)

Create a Fock state for particle number `n ∈ ℕ₀`.
"""
struct Fock <: State
    n::Int64
end

"""
    Crescent(β,ϵ,q)

Create a Cresecent state with parameters `β`, `ϵ`, `q`.
"""
struct Crescent <: State
    β::Complex{Float64}
    ϵ::Complex{Float64}
    q::Float64
end

"""
    Squeezed(β,ϵ)

Create a Squeezed state with parameters `β`, `ϵ`.
"""
struct Squeezed <: State
    β::Complex{Float64}
    ϵ::Complex{Float64}
end

"""
    Squeezed2(r,ϕ)

Create a two-mode squeezed state with parameters `r`, `ϕ`.
"""
struct Squeezed2 <: State
    r::Float64
    ϕ::Float64
end

"""
    Thermal(β,n̄)

Create a Thermal state with parameters `β`, `n̄`.
"""
struct Thermal <: State
    β::Complex{Float64}
    n̄::Float64
end

"""
    Bogoliubov(u,v,n̄)

Create a Bogoliubov state with parameters `u`, `v`, `n̄`.
"""
struct Bogoliubov <: State
    u::Complex{Float64} 
    v::Complex{Float64}
    n̄::Float64
end
