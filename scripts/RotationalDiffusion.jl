using BacterialMotility
using LinearAlgebra
using Distributions
using Plots


#== system geometry ==#
L = 1e3 # box edge (μm)
d = 2 # dimensionality


#== convenience functions ==#
randposition() = (rand(d) .- 1/2) .* (L/2)
function randvelocity(d,U)
    x = rand(d) .- 1/2
    x ./ norm(x) .* U
end # function
randvelocity(U) = randvelocity(d,U)

#== population setup ==#
Δt = 0.1 # s
D_rot = 0.1 # rad²/s
s1 = properties("IntegrationTimestep" => Δt)
s2 = properties("IntegrationTimestep" => Δt,
                "RotationalDiffusivity" => D_rot)
species = ["RT_1", "RT_2"]
N = 30
pop1 = [Bacterium{d}(
    id = species[1], r = randposition(), v = randvelocity(30.0),
    run! = run!, turn! = tumble!, state = copy(s1)) for _ in 1:N]
pop2 = [Bacterium{d}(
    id = species[2], r = randposition(), v = randvelocity(30.0),
    run! = run!, turn! = tumble!, state = copy(s2)) for _ in 1:N]

population = vcat(pop1, pop2)
nbacteria = length(population)
nsteps = 100


#== save routine ==#
function save_position!(traj, bs)
    nbacteria = size(traj, 2)
    t = bs.clock[1]
    for n in 1:nbacteria
        traj[t,n,:] .= bs.population[n].r
    end # for
end # function

trajectories = zeros(nsteps, nbacteria, d)
save_position(bs) = save_position!(trajectories, bs)

#== bacterial system ==#
bs = BacterialSystem(population = population,
                     callback_outer = save_position)

#== simulation ==#
integrate!(bs, nsteps)


#== plot setup ==#
plot_style(palette=:tab10) = (
    thickness_scaling = 1.5,
    guidefontsize = 12,
    tickfontsize = 12,
    legendfontsize = 8,
    grid = false,
    framestyle = :box,
    minorticks = true,
    tick_direction = :in,
    color_palette = palette,
    margin=3Plots.mm,
)

@userplot TailPlot
@recipe function f(p::TailPlot)
    traj, t1, t2, = p.args
    x = @view traj[t1:t2,:,1]
    y = @view traj[t1:t2,:,2]
    seriestype := :path
    linewidth --> 0.5
    label --> false
    x, y
end # recipe

@userplot HeadPlot
@recipe function f(p::HeadPlot)
    traj, t, = p.args
    x = @view traj[t,:,1]
    y = @view traj[t,:,2]
    seriestype := :scatter
    marker --> :circle
    markersize --> 3
    markerstrokewidth --> 0
    label --> false
    x, y
end # recipe

speciescolor = Dict(species .=> 1:length(species))
specieslabels = ["Run-Tumble", "Run-Tumble + Rot. Diffusion"]
linecolor = [speciescolor[bact.id] for _ in 1:1, bact in population]

bgcolor = RGB(0.07, 0.07, 0.07)

#== plot animation ==#
ltail = 35
for t in 2:2:0nsteps
    p = plot(;plot_style(:Dark2)...)
    plot!(p, lims=(-L/2,L/2), aspect_ratio=1, axis=false,
          bgcolor=bgcolor, size=(600,600))
    t0 = max(t-ltail, 1)
    tailplot!(p, trajectories, t0, t; lc=linecolor)
    headplot!(p, trajectories, t; mc=permutedims(linecolor))
    for i in 1:length(species)
        plot!(p, [0.0], [0.0], lab=specieslabels[i],
              lc=speciescolor[species[i]],
              legend=:topright, legendfontsize=5,
              legendfontcolor=:white)
    end # for
    δx = 100 # μm
    xbar0 = -L/2 + L/30
    xbar1 = xbar0 + δx
    ybar = -L/2+L/40
    plot!(p, [xbar0, xbar1], [ybar, ybar], lc=:white, lw=4, lab=false)
    annotate!(p, xbar0+δx/2, ybar+22,
              text("100 μm", :center, :white, 8))
    xtime = L/2 - L/40
    ytime = ybar
    tnow = lpad(round(t*Δt, digits=1), 5, ' ')
    annotate!(p, xtime, ytime,
              text("t = $(tnow) s", :right, :white, 8))
    ndx = lpad(t-1, 4, '0')
    savefig(p, "frame$ndx.png")
end # for
