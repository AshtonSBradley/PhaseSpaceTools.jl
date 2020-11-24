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
