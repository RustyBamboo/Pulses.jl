using Pulses
using Plots; gr()

plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing, grid=true, label="")
scalefontsizes(1.3)

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

sol = Pulses.find_pulse(target, system, Δt, inital_pulse)

# Define loss function
loss(p) = Pulses.compute_loss(target, system, Δt, p)

# 1-Dimensional Linear Interpolation
inter1d_pulse(α) = (1-α) * inital_pulse + α*sol
α_r = LinRange(-50, 50, 1000)
losses = loss.(inter1d_pulse.(α_r))

plot(α_r, losses)

# 2-Dimensional Interpolation
a = rand(100,2)
b = rand(100,2)


inter2d_pulse(α,β) = sol + α*a + β*b
α_r = LinRange(-π, π, 200)
β_r = α_r

P = collect(Iterators.product(α_r, β_r))

losses2d = loss.(inter2d_pulse.(P))

# Plot an animation of the loss landscape
n = 100
animation = @animate for i in range(0, stop = 2π, length = n)
    p = surface(α_r, β_r, losses2d, c=:plasma, showaxis=false, legend=false)
    plot!(p[1], camera = (10 * (5 + 2*cos(i)), 10), background_color = :transparent)
end

gif(animation, "images/fun_loss_landscape.gif", fps=15)
