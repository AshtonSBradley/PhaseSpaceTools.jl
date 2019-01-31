# add real noises for complex variable SDE's
# in DifferentialEquations.jl

#Define real noise for adaptive solver
function realnoise(rand_vec,W,dt,rng)
for i in eachindex(rand_vec)
    rand_vec[i] = randn(rng)
    end
    rand_vec .*= sqrt(abs(dt))
end

function realbridge(rand_vec,W,W0,Wh,q,h,rng)
for i in eachindex(rand_vec)
    rand_vec[i] = randn(rng)
    end
    rand_vec .= sqrt((1 .- q).*q.*abs(h)).*rand_vec.+q.*Wh
end
