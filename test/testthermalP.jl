
nbar = 50
N = 10000
a,ā = thermal(0,nbar,N;dist=:P)
av_a = mean(a);absa = abs(av_a)
n̄ = real(mean(a.*ā))

#TODO tests
