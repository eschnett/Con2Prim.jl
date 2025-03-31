module Con2Prim

using LinearAlgebra
using NonlinearSolve
using StaticArrays

include("solve1d.jl")
include("solve2d.jl")

# internal energy: relativistic ideal gas
const γ = 4 / 3
# See Sackur–Tetrode equation
U(ρ, S, Ne) = ρ^(γ + 1) / (γ + 1) * exp(S / γ) + Ne
dUdρ(ρ, S, Ne) = (γ + 1) / ρ * U(ρ, S, Ne)
dUdS(ρ, S, Ne) = 1 / γ * U(ρ, S, Ne)
dUdNe(ρ, S, Ne) = 1.0

ϵ(ρ, S, Ne) = U(ρ, S, Ne)
press(ρ, S, Ne) = dUdρ(ρ, S, Ne) * ρ^2
temp(ρ, S, Ne) = dUdS(ρ, S, Ne)
μ(ρ, S, Ne) = dUdNe(ρ, S, Ne)

# specific enthalpy
h(ρ, S, Ne) = 1 + ϵ(ρ, S, Ne) + press(ρ, S, Ne) / ρ

# velocity and gamma-factor as function of rapidity w
v(w) = tanh(w)
W(w) = cosh(w)

dens(ρ, S, Ne, w) = ρ * W(w)
dense(ρ, S, Ne, w) = Ne * W(w)
mom(ρ, S, Ne, w) = ρ * h(ρ, S, Ne) * W(w)^2 * v(w)
tau(ρ, S, Ne, w) = ρ * h(ρ, S, Ne) * W(w)^2 - press(ρ, S, Ne) # - dens(ρ, S, Ne, w)

export prim2con
function prim2con(ρ, S, Ne, w)
    return dens(ρ, S, Ne, w), dense(ρ, S, Ne, w), mom(ρ, S, Ne, w), tau(ρ, S, Ne, w)
end
function prim2con(prims::SVector{4})
    ρ, S, Ne, w = prims
    ρ = exp2(ρ)
    # S = exp2(S)
    D, De, M, T = prim2con(ρ, S, Ne, w)
    cons = SVector(D, De, M, T)
    return cons
end

export con2prim
function con2prim(D, De, M, T)
    cons = SVector(D, De, M, T)
    prims = con2prim(cons)
    ρ, S, Ne, w = prims
    ρ = exp2(ρ)
    # S = exp2(S)
    return ρ, S, Ne, w
end
function con2prim(cons::SVector{4})
    ρ, S, Ne, w = (0.9, 0.9, 0.0, 0.0)
    ρ = log2(ρ)
    # S = log2(S)
    prims0 = SVector(ρ, S, Ne, w)
    function f(prims1, extra)
        cons1 = prim2con(prims1)
        return cons1 - cons
    end
    prob = NonlinearProblem(f, prims0, cons)
    solver = SimpleTrustRegion(; autodiff=AutoForwardDiff())
    sol = solve(prob, solver; abstol=1.0e-14, show_trace=Val(true))
    @show sol.u sol.resid sol.retcode sol.stats
    prims = sol.u
    @show prim2con(prims) - cons
    return sol.u
end

end
