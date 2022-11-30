using LinearAlgebra
using Plots; unicodeplots()

function solve(H::Vector{Matrix{ComplexF64}}, Δt::Float64)
    Ω = Δt .* H

    U_i = exp.(-1im * Ω)
    display([norm(u, 2) for u in U_i])

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

display(H[1])
display(H'[1])

display(exp(-1im * H[1]))


U = solve(H, 0.404)


targ = [1. 1.; 1. 1. + 0im]/ 2



# f = norm.([targ - u for u in U], Inf)

# https://journals.aps.org/pra/abstract/10.1103/PhysRevA.72.052337
# f = 2 * 2 .- 2 * real(tr.([targ' * u for u in U]))

# https://arxiv.org/pdf/1101.3817.pdf
# f = real(tr.([targ' * u for u in U])) / 2

# f = tr.([targ' * u / sqrt(2) for u in U])
# f = tr(targ' * U[end])

# https://arxiv.org/pdf/2003.10132.pdf
f = abs2.(tr.([targ' * u for u in U]))

# norms = [norm(u) for u in U]
# display(norms)


# display(U[end])
# display(targ' * U[end])

println(f)

display(plot(x, f))


