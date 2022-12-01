
module BlochSphere
using Makie, CairoMakie, ColorSchemes
set_theme!(theme_black())

function bloch_sphere!(scene)
    wireframe!(scene, Sphere(Point3f(0), 1f0), color=:white, linewidth=0.1, thickness = 0.6f0, transparency = true)
    return scene
end

function bloch_coords(r)
    x = real(r[2,1] + r[1,2])
    y = imag(r[2,1] - r[1,2])
    z = real(r[2,2] - r[1,1])
    return [-x,y,z]
end

function plot_arrow!(scene, r::Matrix{ComplexF64}; kwargs...)
    x, y, z = bloch_coords(r)
    lines!(scene, [0.0, x], [0.0, y], [0.0, z]; kwargs...)
    return scene
end

function plot_points!(scene, rhos::Vector{Matrix{ComplexF64}}; kwargs...)
    out = [bloch_coords(r) for r in rhos]
    x = getindex.(out, 1)
    y = getindex.(out, 2)
    z = getindex.(out, 3)
    scatter!(scene, x, y, z, color=1:length(out), colormap=:plasma, markersize=20)
    # lines!(scene, [0.0, x], [0.0, y], [0.0, z]; kwargs...)
    return scene
    # x,y,z
end


function create_bloch()
    scene = Scene(resolution=(300, 300), center=false)
    # Makie.scale!(scene, -3,-3,-3)
    Makie.rotate!(scene, Vec3f(1,0,0), pi/3)
    bloch_sphere!(scene)
    scene
end


function save(filename, scene)
    Makie.save(filename, scene)
end

# plot_arrow!(scene, [1. 0; 0 0im])
end