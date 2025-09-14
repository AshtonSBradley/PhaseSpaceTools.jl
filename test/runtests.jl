#!/usr/bin/env julia

using PhaseSpaceTools, Test

function randuv()
    u,v = randnc(2)
    nrm = abs2(u)-abs2(v)
    if nrm < 0
        new = [0 1;1 0]*[u; v]
        u,v = new
    end
    u /=sqrt(abs(nrm))
    v /=sqrt(abs(nrm))
    return u,v
end

@testset "Arg tests" begin 

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

end

@testset "Noises" begin 

    N = 100_000
    a = randnc(N)
    mean(abs2.(a))

    @test isapprox(mean(abs2.(a)),1.0,rtol=5e-2)
    @test isapprox(abs(mean(a)),0.0,atol=1e-2)
end

@testset "Coherent W" begin 

    # Coherent
    N = 100000
    α = 100
    state = Coherent(α)
    a,ā = wigner(state,N)

    meana = mean(a)
    n̄ = real(mean(@. ā*a)-.5)
    Vn = mean(@. ā^2*a^2)-mean(@. a*ā)-n̄^2 |> abs

    #test
    @test norm(meana - α)/abs(α) < 0.01
    @test abs(n̄ - abs(α).^2)/abs(α)^2 < 0.01
    @test abs(Vn - abs(α)^2)/abs(α)^2 < 0.05

end

@testset "Coherent +P" begin 

    N = 1000
    α = 100*randnc()
    state = Coherent(α)
    a,ā = positiveP(state,N)


    meana = mean(a)
    n̄ = real(mean(@. ā*a))
    Vn = mean(@. ā^2*a^2)-mean(@. a*ā)-n̄^2 |> abs

    #test
    @test n̄ ≈ abs2(α)
    @test norm(meana - α)/abs(α) < 0.01
    @test abs(n̄ - abs(α).^2)/abs(α)^2 < 0.01
    @test abs(Vn - abs(α)^2)/abs(α)^2 < 0.05

end


@testset "Thermal P,Q,W" begin 

    N = 100001
    n = 50
    β = 0.
    s = Thermal(β,n)

    a,ā = glauberP(s,N)
    na = real(mean(a.*ā))
    @test isapprox(na,n,atol=.9)

    a,ā = wigner(s,N)
    na = real(mean(a.*ā))-0.5
    @test isapprox(na,n,atol=.9)

    a,ā = husimiQ(s,N)
    na = real(mean(a.*ā))-1.0
    @test isapprox(na,n,atol=.9)

end

@testset "Fock W" begin 

    # Fock
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

end

@testset "Fock +W" begin 

#test +W

    n = 99
    N = 100_000
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

end

@testset "Fock +W (n≫1)" begin 

#test +W for large n

    n = 101
    state = Fock(n)
    N = 100000

    a,ā = positiveW(state,N)

    av_a = mean(a)
    absa = abs(av_a)
    n̄ = real(mean(a.*ā)-.5)
    Vn= abs(mean(a.*a.*ā.*ā)-mean(a.*ā)-n̄.^2)
    rel_num_var = sqrt(abs(Vn))/abs(n̄)

    #test
    @test absa < 0.1
    @test (n̄ - n)/n < 0.01
    @test Vn < 30
    @test rel_num_var < 0.1

end

@testset "Fock +P" begin 

    n = 99
    N = 100_000
    state = Fock(n)
    a,ā = positiveP(state,N)

    meana = mean(a)
    absa = abs(meana)
    n̄ = mean(a.*ā) |> real
    Vn = abs(mean(a.*a.*ā.*ā))-abs(mean(a.*ā)).^2
    rel_num_var = sqrt(abs(Vn))/abs(n̄);

    #test
    @test absa < 0.1
    @test n̄ - n < 1
    @test Vn < 50
    @test rel_num_var < 0.15

end

@testset "Squeezed +P " begin 

    #Squeezed
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
    @test norm(av_a - β)/abs(β) < 0.02
    @test abs(n̄ - nbar)/abs(β)^2 < 0.02

end

# TODO:
@testset "Squeezed2 +P " begin 

## Squeezed
r = .5
ϕ = 0
n̄ = sinh(r)^2
state = SqueezedTwoMode(r,ϕ)
N = 1000000

a,a⁺,b,b⁺ = positiveP(state,N)

na = mean(a.*a⁺) |> real
nb = mean(b.*b⁺) |> real
@test isapprox(na,n̄,atol=5e-2)
@test isapprox(nb,n̄,atol=5e-2) 

##quadratures
X = a + b⁺
X⁺ = a⁺ + b
Y = im*(a⁺ - b)
Y⁺ = -im*(a - b⁺)
σX = mean(X.*X⁺)+1 |> real |> sqrt
σY = mean(Y.*Y⁺)+1 |> real |> sqrt

@test isapprox(σX,exp(r),rtol=5e-2)
@test isapprox(σY,exp(-r),rtol=5e-2)

##TODO ϕ=π/4 test
end

@testset "Squeezed W" begin 

    # Squeezed
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

end

@testset "Bogoliubov W" begin 

    # Bogoliubov
    N = 100000
    n̄ = 10

    u,v = randuv()
    abs2(u)-abs2(v) ≈ 1.0

    #thermal state
    state = Bogoliubov(u,v,n̄)
    a,ā = wigner(state,N)

    # particle mode population
    na = real(mean(a.*ā)) - 0.5
    nth = (abs2(u) + abs2(v))*(n̄ + 0.5) - 0.5
    @test isapprox(na,nth,rtol=1e-1)

    #test vacuum limit
    state = Bogoliubov(u,v,0)
    a,ā = wigner(state,N)

    na = real(mean(a.*ā)) - 0.5
    nvac = abs2(v)
    @test isapprox(na,nvac,atol=0.3)

end

@testset "Bogoliubov Q" begin 

    # Bogoliubov
    N = 100000
    n̄ = 10

    u,v = randuv()
    abs2(u)-abs2(v) ≈ 1.0

    #thermal state
    state = Bogoliubov(u,v,n̄)
    a,ā = husimiQ(state,N)

    # particle mode population
    na = real(mean(a.*ā)) - 1.0
    nth = (abs2(u) + abs2(v))*(n̄ + 1.0) - 1.0
    @test isapprox(na,nth,rtol=1e-1)

    #test vacuum limit
    state = Bogoliubov(u,v,0)
    a,ā = husimiQ(state,N)

    na = real(mean(a.*ā)) - 1.0
    nvac = abs2(v)
    @test isapprox(na,nvac,atol=0.3)

end

@testset "Bogoliubov +P" begin 

    # Bogoliubov
    N = 100000
    n̄ = 10

    u,v = randuv()
    @test abs2(u)-abs2(v) ≈ 1.0

    #thermal state
    state = Bogoliubov(u,v,n̄)
    a,ā = positiveP(state,N)

    # particle mode population
    na = mean(a.*ā) |> real
    nth = (abs2(u) + abs2(v))*n̄ + abs(v)
    @test isapprox(na,nth,rtol=1e-1)

    #test vacuum limit
    state = Bogoliubov(u,v,0)
    a,ā = positiveP(state,N)

    na = mean(a.*ā) |> real
    nvac = abs2(v)
    @test isapprox(na,nvac,atol=0.3)

end


