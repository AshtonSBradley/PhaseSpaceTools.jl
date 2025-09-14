using CairoMakie, LaTeXStrings
using ClassicalOrthogonalPolynomials
# set_theme!(theme_black())
set_theme!(theme_latexfonts())

"""
    phi_N = hermite_mode(x,N)

Physicists Hermite mode with Gaussian weight at `x`. Normalisation ∫|ψ|²dx = 1. 
Uses stable recursion.
"""
function hermite_mode(x,N)
    Z₀ = exp(-x^2/2)/π^(0.25)
    Z₁ = √2*x*Z₀
    if N == 0 
        return Z₀
    elseif N == 1
        return Z₁
    else
        Zₙ₊₁,Zₙ,Zₙ₋₁ = zero(x),Z₁,Z₀
        for n in 1:N-1
            Zₙ₊₁ = sqrt(2/(n+1))*x*Zₙ-sqrt(n/(n+1))*Zₙ₋₁
            Zₙ,Zₙ₋₁ = Zₙ₊₁,Zₙ
        end
    return Zₙ₊₁
    end
end

W(x,p,n) = (-1)^n/π * exp(-(x^2+p^2)) * laguerrel(n, 2*(x^2 + p^2))
ψ(x,n) = hermite_mode(x,n)
P(x,n) = abs2(ψ(x,n))

## Plot it!
n = 15
Nx = 1000
xmax = sqrt(2n+2)+2
x = range(-xmax,xmax,Nx)
p = x

f = Figure(size=(600,600))
# cmap = reverse ? Reverse(cmap) : cmap
s = [W(x,p,n) for x in x, p in p]

ax = Axis3(f[1,1]; aspect =(1,1,1),
perspectiveness = .6f0,
azimuth = -1.275π * 1.81,
elevation = pi/6, protrusions=0,
xlabel=L"x", ylabel=L"p", 
xlabelalign=(:right,:center),
ylabelalign=(:left,:center),
xlabelrotation=0,ylabelrotation=0,
zlabel="",
xticks=LinearTicks(3),yticks=LinearTicks(3),zticks=LinearTicks(3),
)
hidedecorations!(ax)
surface!(ax,x,p,s,colormap=:cool)
f


##
save("fock.pdf", f)