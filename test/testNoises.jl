N = 100_000
a = randnc(N)
mean(abs2.(a))

@test isapprox(mean(abs2.(a)),1.0,rtol=5e-2)
