using Pulses
using Plots; gr()
using LaTeXStrings
include("../src/bloch_sphere.jl")
using Main.BlochSphere
theme(:dracula)

# Drift Hamiltonian
H_0 = [0. 0. 0.;
       0. 0. 0.;
       0. 0. -1.]

# Control Hamiltonians       
X = [0. 1. 0.;
     1. 0. sqrt(2);
     0 sqrt(2) 0]

Y = [0. 1im 0.;
     -1im 0. 1im*sqrt(2);
      0 -1im*sqrt(2) 0]
      
# Desired Unitary operator
target = [1. 1. 0;
          1. -1 0;
          0im 0 sqrt(2)
          ] / sqrt(2)

# Define our system using drift and control Hamiltonians
system = Pulses.System(H_0, [X, Y])

# Define a time scale
t_r = LinRange(0, 10, 100)
Δt = t_r[2] - t_r[1]

# Define an initial pulse
x = LinRange(-3, 3, 100)
inital_pulse = 1e-3.*[exp.(-x.^2) exp.(-x.^2)]

# Solve for an optimal pulse
sol = Pulses.find_pulse(target, system, Δt, inital_pulse)


# ------------ Plotting -------------- #

# Plot the pulse
plot(layout=(2,1))
plot!(subplot=1, t_r, sol, label=["x" "y"])
plot!(subplot=1, title="Pulses", xlabel="Time [ns]", ylabel="Amplitude")

# Plot the state evolution with optimal pulse
H = Pulses.hamiltonian(system, sol)
evolution = Pulses.solve_state_history(H, Δt, [1; 0im; 0.])

p_0 = [abs2(e[1]) for e in evolution]
p_1 = [abs2(e[2]) for e in evolution]

plot!(subplot=2, t_r, p_0, label=L"$|0\rangle$")
plot!(subplot=2, t_r, p_1, label=L"$|1\rangle$")
plot!(subplot=2, xlabel="Time [ns]", ylabel="Probability")

savefig("images/hadamard.svg")

# Plot on the Bloch sphere

rhos = [e[begin:2] * e[begin:2]' for e in evolution]

s = BlochSphere.create_bloch()
s = BlochSphere.plot_arrow!(s, [1. 0; 0 0im], color=:green)
s = BlochSphere.plot_arrow!(s, [1. 1; 1 1+0im]/2, color=:red)
s = BlochSphere.plot_points!(s, rhos)
# display(s.scene)
BlochSphere.save("images/hadamard-bloch.png", s)