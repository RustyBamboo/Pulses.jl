using LinearAlgebra
using Zygote
using Plots; unicodeplots()

function solve(H::Vector{Matrix{ComplexF64}}, Δt::Float64)
    Ω = Δt .* H

    U_i = exp.(-1im * Ω)

    U = cumprod(U_i)

    return U
end


function solve_state(H::Vector{Matrix{ComplexF64}}, Δt::Float64, ψ0::Vector{ComplexF64})
    U = I
    Ω = Δt .* H

    U_i = exp.(-1im * Ω)

    U = cumprod(U_i)

    ψ_t = [u * ψ0 for u in U] 

    return ψ_t
end

H_0 = [0. 0.; 0. 0.]
X = [0. 1.; 1. 0.]
Y = [0. -1im; 1im 0.]

t_r = LinRange(0, 40, 100)
x = LinRange(-3, 3, 100)
pulse = exp.(-x.^2)

H = [p * X for p in pulse] + [p * Y for p in pulse]
U = solve(H, 0.404)


function compute(params)
    H = [p * X for p in params] + [p * Y for p in params]
    U = solve(H, 0.404)
    abs2.(tr.([targ' * u for u in U]))
end

j = jacobian(a -> compute(a), pulse)

# TODO: fix Hessian
# h = hessian(a -> compute(a), pulse)


