using PhaseSpaceTools, Test, Statistics

## Two-mode Squeezed states
r = 1
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

