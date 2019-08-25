
# Coherent
state = Coherent(10.0)
@test typeof(state) <: State
@test typeof(wigner(state,1)[1][1]) == Complex{Float64}
@test typeof(glauberP(state,1)[1][1]) == Complex{Float64}
@test typeof(husimiQ(state,1)[1][1]) == Complex{Float64}
@test typeof(positiveP(state,1)[1][1]) == Complex{Float64}
