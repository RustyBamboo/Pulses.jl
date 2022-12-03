module BlochSphere
using Makie, CairoMakie, ColorSchemes, Colors
set_theme!(theme_black())



spherecolor = function bloch_sphere!(scene,
    spherecolor=RGBA(0.5, 0.5, 0.5, 0.8),
    wirecolor=RGBA(1.0, 1.0, 1.0, 0.2),
    wirecolor_outline=RGBA(1.0, 1.0, 1.0, 0.8))

    r = 1.0
    function sphere(nu, mu, u1, u2, r=1.0)
        u = range(u1, u2; length=nu)
        v = range(0, π; length=mu)
        x = r * cos.(u) * sin.(v)'
        y = r * sin.(u) * sin.(v)'
        z = r * ones(nu) * cos.(v)'
        return (x, y, z)
    end

    surface!(scene, sphere(25, 25, 0, 2 * π, r)..., color=fill(spherecolor, 25, 25), transparency=true)
    # surface!(scene, sphere(25, 25, -π, 0)..., color=fill(RGBA(1.,1.,1.,0.8),25,25), transparency = true)
    wireframe!(scene, sphere(13, 13, 0, 2 * π, r)..., linewidth=1, color=wirecolor, transparency=true)

    n = 25
    u = range(0, 2 * π; length=n)
    v = range(0, π; length=n)
    lines!(scene, r * cos.(u), r * sin.(u), r * zeros(n), color=wirecolor_outline)
    lines!(scene, r * zeros(n), r * sin.(u), r * cos.(u), color=wirecolor_outline)
    lines!(scene, r * sin.(u), r * zeros(n), r * cos.(u), color=wirecolor_outline)
    lines!(scene, [0.0, 0.0, 0], [0.0, r, 0], [0, 0, 0.0], color=wirecolor_outline)
    lines!(scene, [0.0, 0.0, 0], [0.0, -r, 0], [0, 0, 0.0], color=wirecolor_outline)
    lines!(scene, [0.0, 0.0, 0], [0.0, 0.0, 0], [0, 0, r], color=wirecolor_outline)
    lines!(scene, [0.0, 0.0, 0], [0.0, 0.0, 0], [0, 0, -r], color=wirecolor_outline)
    lines!(scene, [0.0, 0.0, r], [0.0, 0.0, 0], [0, 0, 0.0], color=wirecolor_outline)
    lines!(scene, [0.0, 0.0, -r], [0.0, 0.0, 0], [0, 0, 0.0], color=wirecolor_outline)

    return scene
end


function bloch_coords(r)
    x = real(r[2, 1] + r[1, 2])
    y = imag(r[2, 1] - r[1, 2])
    z = real(r[2, 2] - r[1, 1])
    return [x, y, -z]
end

function plot_arrow!(scene, r::Matrix{ComplexF64}; kwargs...)
    x, y, z = bloch_coords(r)
    # lines!(scene, [0.0, x], [0.0, y], [0.0, z]; kwargs...)
    arrows!([0.0], [0.0], [0.0], [x], [y], [z], linewidth=0.03f0, arrowsize=0.1f0, lengthscale=0.9; kwargs...)
    return scene
end

function plot_points!(scene, rhos::Vector{Matrix{ComplexF64}}; kwargs...)
    out = [bloch_coords(r) for r in rhos]
    x = getindex.(out, 1)
    y = getindex.(out, 2)
    z = getindex.(out, 3)
    scatter!(scene, x, y, z, color=1:length(out), colormap=:plasma, markersize=10)
    # lines!(scene, [0.0, x], [0.0, y], [0.0, z]; kwargs...)
    return scene
end


function create_bloch()
    fig = Figure(resolution=(500, 500), backgroundcolor=RGBf(0.7, 0.8, 1))
    scene = LScene(fig[1, 1], show_axis=false,)
    bloch_sphere!(scene)
    return scene
end

function save(filename, scene)
    return Makie.save(filename, scene.scene)
end

end
