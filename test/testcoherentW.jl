#test W

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
