"""
    x = reject(P,w,N,Pmax)

Generate `x` distributed according to a probability distribution by rejection sampling over a finite window.

`P`: probability distribution, normalized to 1.

`w=[w1,w2]`: window for sampling `x`.

`N`: number of samples.

`Pmax`: numerical upper bound for `P`: `Pmax ≧ max(P(w))`.

"""
function reject(P,w,N,Pmax)
    samples = Array{Float64}(undef,N)
    count = 1
    while count < N+1
        y = w[1] + rand()*(w[2] - w[1])
        z = rand()*Pmax
        (z < P(y)) && (samples[count] = y)
        count += 1
    end
    return samples
end
