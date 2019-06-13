abstract type State end

struct Coherent <: State
    β::Complex{Float64}
end

struct Fock <: State
    n::Int64
end

struct Crescent <: State
    β::Complex{Float64}
    ϵ::Complex{Float64}
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
