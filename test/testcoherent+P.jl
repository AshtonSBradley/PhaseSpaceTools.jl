N = 1000
α = 100*crandn()
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
