using LinearAlgebra
using Zygote
using Optim
using Plots; gr()

function solve(H::Vector{Matrix{ComplexF64}}, Δt::Float64)
    Ω = Δt .* H

    U_i = exp.(-1im * Ω)

    U = prod(U_i)

    return U
end


function solve_state_history(H::Vector{Matrix{ComplexF64}}, Δt::Float64, ψ0::Vector{ComplexF64})
    U = I
    Ω = Δt .* H

    U_i = exp.(-1im * Ω)

    U = cumprod(U_i)

    ψ_t = [u * ψ0 for u in U] 

    return ψ_t
end

H_0 = [0. 0. 0.; 0. 0. 0.; 0. 0. -1.]
X = [0. 1. 0.; 1. 0. sqrt(2); 0 sqrt(2) 0]
# Y = [0. -1im; 1im 0.]
Y = [0. 1im 0.; -1im 0. 1im*sqrt(2); 0 -1im*sqrt(2) 0]


t_r = LinRange(0, 40, 100)
x = LinRange(-3, 3, 100)
pulse = 1e-3.*[exp.(-x.^2) exp.(-x.^2)]

targ = [1. 1. 0; 1. -1+0im 0; 0 0 sqrt(2)] / sqrt(2)

H = [H_0 + p * X for p in pulse[:,1]] + [p * Y for p in pulse[:,2]]

X
Y

function compute(params)
    H = [H_0 + p * X for p in params[:,1]] + [p * Y for p in params[:,2]]
    U = solve(H, 0.404)
    # abs2.(tr.([targ' * u for u in U]))
    # abs2.(tr(targ' * U[end]))
    # norm(targ - U[end])^2
    u_targ_norm_psu = abs2(tr(targ * U))
    1 - u_targ_norm_psu
end

function compute_j(params)
    # jacobian(a -> compute(a), params)
    reshape(jacobian(a -> compute(a), params)[1], (100, 2))

end

function compute_u(params)
    H = [H_0 + p * X for p in params[:,1]] + [p * Y for p in params[:,2]]
    solve(H, 0.404)
end

# result = optimize(compute,  pulse, LBFGS())
# result = optimize(compute,  pulse, LBFGS())

result = optimize(compute, compute_j, pulse, LBFGS(); inplace = false)

sol = Optim.minimizer(result)

# plot(x, pulse)
plot(t_r, sol, label=["x" "y"])
plot!(title="Pulses", xlabel="Time [ns]", ylabel="Amplitude")
savefig("hadamard.png")

Optim.minimum(result)

targ

U_f = compute_u(sol)