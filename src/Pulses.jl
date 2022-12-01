module Pulses

using LinearAlgebra
using Zygote
using Optim

"""
    System(drift, control)

Stores information about the drift Hamiltonian and all the control Hamiltonians
"""
struct System
    drift::Matrix{ComplexF64}
    control::Vector{Matrix{ComplexF64}}
end

"""
    hamiltonian(s, pulses)

Takes a System (drift+control Hamiltonians) and pulses to create a time-dependent Hamiltonian
"""
function hamiltonian(s::System, pulses::Matrix{Float64})
    @assert size(s.control, 1) == size(pulses, 2) "Control parameters and size of pulses do not match"
    # Zygote did not like the adjoint:
    # H = Ref(s.drift) .+ sum(s.control' .* pulses, dims = 2)
    # so instead manually do the shape matching for broadcasting
    H = Ref(s.drift) .+ sum([Ref(s.control[i]) .* pulses[:, i] for i = 1:size(pulses, 2)])
    return H
end

"""
    solve(H, Δt)

First-order solution to the time-dependent Schrodinger equation
"""
function solve(H::Vector{Matrix{ComplexF64}}, Δt::Float64)
    Ω = Δt .* H
    U_i = exp.(-1im * Ω)
    U = prod(U_i)
    return U
end

"""
    solve_state_history(H, Δt, ψ0)

First-order solution to the time-dependent Schrodinger equation, but keeps a history of how a state ψ0 evolves in time
"""
function solve_state_history(H::Vector{Matrix{ComplexF64}}, Δt::Float64, ψ0::Vector{ComplexF64})
    U = I
    Ω = Δt .* H
    U_i = exp.(-1im * Ω)
    U = cumprod(U_i)
    ψ_t = [u * ψ0 for u in U]
    return ψ_t
end

"""
    Solve the system and compute loss with respect to a target
"""
function compute_loss(target::Matrix{ComplexF64}, s::System, Δt::Float64, params)
    H = hamiltonian(s, params)
    U = solve(H, Δt)
    u_targ_norm_psu = abs2(tr(target' * U))
    1 - u_targ_norm_psu
end

"""
    Compute the gradient of the loss
"""
function compute_j_loss(target::Matrix{ComplexF64}, s::System, Δt::Float64, params)
    reshape(jacobian((a) -> compute_loss(target, s, Δt, a), params)[1], size(params))
end

"""
    Get final U
"""
function compute_u(s::System, pulse, Δt::Float64)
    H = hamiltonian(s, pulse)
    return solve(H, Δt)
end

function find_pulse(target::Matrix{ComplexF64}, s::System, Δt::Float64, initial_pulse)
    loss(x) = compute_loss(target, s, Δt, x)
    j_loss(x) = compute_j_loss(target, s, Δt, x)
    result = optimize(loss, j_loss, initial_pulse, LBFGS(); inplace = false)
    sol = Optim.minimizer(result)
    sol
end




