"""
    x = reject(P,w,N,Pmax)

Generate `x` distributed according to a probability distribution by rejection sampling over a finite window.

`P`: probability distribution, normalized to 1.

`w=[w1,w2]`: window for sampling `x`.

`N`: number of samples.

`Pmax`: numerical upper bound for `P`: `Pmax ≧ max(P(w))`.

"""
function reject(P,w,N,Pmax)
    samples = Array{Float64}(1)
    while length(samples) < N+1
        y = w[1] + rand()*(w[2]-w[1])
        z = rand()*Pmax
        z < P(y) && push!(samples,y)
    end
    return samples[2:end]
end
