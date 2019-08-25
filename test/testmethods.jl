# Coherent
state = Coherent(10.0)
@test typeof(state) <: State
@test typeof(wigner(state,1)[1][1]) == Complex{Float64}
@test typeof(glauberP(state,1)[1][1]) == Complex{Float64}
@test typeof(husimiQ(state,1)[1][1]) == Complex{Float64}
@test typeof(positiveP(state,1)[1][1]) == Complex{Float64}

# Fock
n = 100
state = Fock(n)
@test typeof(state) == typeof(Fock(n))
@test typeof(state) <: State
@test typeof(wigner(state,1)[1][1]) == Complex{Float64}
@test typeof(husimiQ(state,1)[1][1]) == Complex{Float64}
@test typeof(positiveP(state,1)[1][1]) == Complex{Float64}

# Bogoliubov
state = Bogoliubov(.1,.1,1)
@test typeof(state) <: State
@test typeof(state) == typeof(Bogoliubov(.1,.1,1))
@test typeof(wigner(state,1)[1][1]) == Complex{Float64}
@test typeof(glauberP(state,1)[1][1]) == Complex{Float64}
@test typeof(husimiQ(state,1)[1][1]) == Complex{Float64}
@test typeof(positiveP(state,1)[1][1]) == Complex{Float64}

# Squeezed
state = Squeezed(randnc(),randnc())
@test typeof(state) <: State
@test typeof(state) == typeof(Squeezed(randnc(),randnc()))
@test typeof(wigner(state,1)[1][1]) == Complex{Float64}
@test typeof(husimiQ(state,1)[1][1]) == Complex{Float64}
@test typeof(positiveP(state,1)[1][1]) == Complex{Float64}

# Thermal
state = Thermal(.1,1)
@test typeof(state) <: State
@test typeof(state) == typeof(Thermal(.1,1))
@test typeof(wigner(state,1)[1][1]) == Complex{Float64}
@test typeof(glauberP(state,1)[1][1]) == Complex{Float64}
@test typeof(husimiQ(state,1)[1][1]) == Complex{Float64}
@test typeof(positiveP(state,1)[1][1]) == Complex{Float64}

# Crescent
state = Crescent(1.0,0.1,0.1)
@test typeof(state) <: State
@test typeof(state) == typeof(Crescent(1.0,0.1,0.1))
@test typeof(wigner(state,1)[1][1]) == Complex{Float64}
@test typeof(husimiQ(state,1)[1][1]) == Complex{Float64}
@test typeof(positiveP(state,1)[1][1]) == Complex{Float64}

# Helpers
f(x) = exp(-x^2)
@test typeof(reject(f,[-1,1],1,1.1)[1]) == Float64
@test typeof(plaguerre(.1,1)) == Float64
@test typeof(plaguerre_asymptotic(.1,1)) == Float64
@test typeof(laguerren(.1,1)) == Float64
@test typeof(plaguerre_asymptotic(.1,1)) == Float64
