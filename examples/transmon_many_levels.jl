using Pulses
using Plots;
gr();
using LaTeXStrings
using LinearAlgebra

plot_font = "Computer Modern"
default(fontfamily=plot_font,
    linewidth=2, framestyle=:box, label=nothing, grid=true)

transitions = [4.80637, 4.58210259, 4.33777661]
levels = size(transitions, 1) + 1

# Basic approximation
E_c = (transitions[1] - transitions[2])
E_j = (transitions[2] + E_c)^2 / (8 * E_c)
ng = 0

n_cut = 40

n = diagm(-n_cut:n_cut)

hc_n2 = 4 * n * n
hc_n1 = -2 * 4 * n
hc_n0 = 4 * I
hj = -0.5 * (diagm(1 => ones(2 * n_cut)) + diagm(-1 => ones(2 * n_cut)))

H_T = E_c * (hc_n2 + ng * hc_n1 + ng^2 * hc_n0) + E_j * hj

F = eigen(H_T)

eigs = F.values
evec = F.vectors

N_operator = diagm(-n_cut:n_cut)

charge(i, j) = evec[:, i]' * N_operator * evec[:, j]

charge_hat = reshape([charge(i, j) for i in 1:levels for j in 1:levels], (levels, levels))
charge_hat = charge_hat / charge_hat[1, 2]

c_destroy = triu(charge_hat)
c_create = tril(charge_hat)

tilde_c_dag = diagm(-1 => diag(charge_hat, -1))
tilde_c = diagm(1 => diag(charge_hat, 1))

h_c_0 = tilde_c + tilde_c_dag
h_c_1 = 1im * (tilde_c - tilde_c_dag)

freq = [0; transitions]
wb = cumsum(freq) * 2π

ad_a = diagm(1:levels)
h_d = diagm(wb) - ad_a * wb[2]


system = Pulses.System(h_d, [h_c_0, h_c_1])
t_r = LinRange(0, 10, 100)
Δt = t_r[2] - t_r[1]

x = LinRange(-3, 3, 100)
inital_pulse = 1e-3 .* [exp.(-x .^ 2) exp.(-x .^ 2)]

# desired = [1. 1.; 1. -1.] / sqrt(2)
# target = diagm([1.,1,1,1+0im])
# target[begin:2,begin:2] = desired

target = [1.0 0.0 0.0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0im]

target

sol = Pulses.find_pulse(target, system, Δt, inital_pulse)

plot(layout=(2, 1))
plot!(subplot=1, t_r, sol, label=["x" "y"])
plot!(subplot=1, title="Pulses", xlabel="Time [ns]", ylabel="Amplitude")

final_u = Pulses.compute_u(system, sol, Δt)
abs2.(final_u)
# Landscape

# Define loss function
loss(p) = Pulses.compute_loss(target, system, Δt, p)

function get_landscape()
    a = rand(100, 2)
    b = rand(100, 2)


    inter2d_pulse((α, β)) = inital_pulse + α * a + β * b
    α_r = LinRange(-π, π, 200)
    β_r = α_r

    P = collect(Iterators.product(α_r, β_r))
    losses2d = loss.(inter2d_pulse.(P))
    return losses2d
end

sample_n = 10

losses_all = zeros(200, 200, sample_n)

Threads.@threads for i = 1:sample_n
    losses_all[:, :, i] = get_landscape()
end

using JLD
save("data/transmon_many_levels_rand_landscape.jld", "data", losses_all)
losses_all = load("data/transmon_many_levels_rand_landscape.jld")["data"]

using Statistics
losses_all_mean = Statistics.mean(losses_all, dims=3)

α_r = LinRange(-π, π, 200)
β_r = α_r
surface(α_r, β_r, losses_all_mean[:, :, 1], axis=(yticklabelsvisible=false,), xlabel=L"\alpha_1", ylabel=L"\alpha_2", zlabel=L"\mathcal{F}", legend=false, camera=(60, 10))

savefig("images/cnot_loss.pdf")
# plot!(camera = (10 * (5 + 3*cos(i)), 10))

animation = @animate for i in range(0, stop=2π, length=100)

    p = surface(α_r, β_r, losses_all_mean[:, :, 1],  xlabel=L"\alpha_1", ylabel=L"\alpha_2", zlabel=L"\mathcal{F}", legend=false)
    # plot!(p, title="Loss Landscape")
    plot!(p, camera=(10 * (5 + 3 * cos(i)), 10))
    # plot!(p, showaxis=false)
    plot(p, background_color = :transparent)

end

gif(animation, "images/cnot_plot_loss.gif", fps=15)
webm(animation, "images/cnot_plot_loss.webm", fps=15)

theme(:dark)
using Plots.PlotMeasures
animation = @animate for i in 1:130
    l = @layout [
        a{0.7w} b
    ]

    if i > 100
        idx = 100
    else
        idx = i
    end



    p1 = plot(t_r, sol, label=[L"\alpha_0" L"\alpha_1"])
    plot!(p1, title="Pulses", xlabel="Time [ns]", ylabel="Amplitude")
    vline!(p1, [t_r[idx]], line=(:black, 1), label=:none)
    vspan!(p1, [0, t_r[idx]]; color=:black, alpha=:0.2, label=:none)
    # fixed_y = zeros(40)
    # plot!(p[2], i, )

    final_u = Pulses.compute_u(system, sol[begin:idx, :], Δt)
    p2 = heatmap((abs2.(final_u)),
        yflip=true, colorbar=false, title="Quantum Gate", vmin=0.0, vmax=1.0
        )
    plot!(p2, xticks=([1:4;], [L"|00\rangle", L"|01\rangle", L"|10\rangle", L"|11\rangle"]))
    plot!(p2, yticks=([1:4;], [L"\langle 00|", L"\langle 01|", L"\langle 10|", L"\langle 11|"]))

    plot(p1, p2, layout = l, size=(700, 200), bottom_margin=5mm, left_margin=5mm)
    # plot!(p[3], showaxis=false)
end

# gif(animation, "images/cnot_plot.gif", fps=15)
gif(animation, "images/cnot_plot_gate.gif", fps=15)

webm(animation, "images/cnot_plot_gate.webm", fps=15)


