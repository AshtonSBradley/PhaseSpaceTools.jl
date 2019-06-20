#define all

using Revise, PhaseSpaceTools

import PhaseSpaceTools
instanceof(x::Type) = x.instance

dist = subtypes(PhaseSpace) .|> Symbol

#methods for +W incomplete
dist = [dist[i] for i in (1,2,3,5)]

meth = subtypes(PhaseSpace) .|> instanceof
stat = subtypes(State)

methodname(T) where T<:PhaseSpace
    


coherent(β,N,T) where T = meth[1](Coherent(β),N)


function phasespace(T::PhaseSpace) where T
    return

for (me, di) in zip(meth,dist)
    @eval coherent(β,N,$di) = $me(Coherent(β),N)
end

coherent(β,N,dist::PhaseSpace) = $me(Coherent(β),N)
