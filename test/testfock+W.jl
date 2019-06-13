#test +W

n = 99
N = 100000
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
