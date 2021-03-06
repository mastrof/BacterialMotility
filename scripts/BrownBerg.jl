using BacterialMotility
using LinearAlgebra
using Distributions
using Plots


#== system geometry ==#
L = 1000 # μm
d = 3 # dimensionality

#== convenience functions ==#
randposition() = (rand(d) .- 1/2) .* (L/2)
function randvelocity(d,U)
    x = rand(d) .- 1/2
    x ./ norm(x) .* U
end # function
randvelocity(U) = randvelocity(d,U)


#== motile patterns ==#
custom_revflick!(bacterium) = reverse_flick!(bacterium, Normal(π,π/8), Normal(π/2,π/8))


#== field setup ==#
R = 20.0
C = 5.0
Cbg = 0.01
f = Concentration_SteadyDiffusionSphericalSource_3D(R=R, C=C, Cbg=Cbg)


#== population setup ==#
Δt = 0.067 # s
D_rot = 0.0325 # rad²/s
s1 = propertiesBrownBerg("IntegrationTimestep" => Δt,
                         "RotationalDiffusivity" => D_rot,
                         "MotorGain" => 0.0)
s2 = propertiesBrownBerg("IntegrationTimestep" => Δt,
                        "RotationalDiffusivity" => D_rot)

species = ["Random Walk", "Chemotaxis"]
N = 50
pop1 = [BacteriumBrownBerg{d}(
    id = species[1], r = randposition(), v = randvelocity(36.0),
    state = copy(s1)) for _ in 1:N]
pop2 = [BacteriumBrownBerg{d}(
    id = species[2], r = randposition(), v = randvelocity(36.0),
    state = copy(s2)) for _ in 1:N]

population = vcat(pop1, pop2)
nbacteria = length(population)
nsteps = 1000

#== callbacks ==#
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
                     field = f,
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
specieslabels = species
linecolor = [speciescolor[bact.id] for _ in 1:1, bact in population]

function cfield(x,y)
    r = sqrt(x*x+y*y)
    r < R ? Cbg+C : Cbg+C*R/r
end # function
clims = (cfield(L/2,L/2), cfield(0,0))
xx = yy = -L/2:4:L/2
cc = cfield.(xx',yy)

bgcolor = RGB(0.07, 0.07, 0.07)

#== plot animation ==#
ltail = 25
for t in 2:4:0nsteps
    p = plot(;plot_style(:Dark2)...)
    plot!(p, lims=(-L/2,L/2), aspect_ratio=1, axis=false,
          bgcolor=bgcolor,  size=(600,600))
    t0 = max(t-ltail, 1)
    heatmap!(p, xx, yy, cc, c=:bone, clims=clims, cbar=false)
    tailplot!(p, trajectories, t0, t; lc=linecolor)
    headplot!(p, trajectories, t; mc=permutedims(linecolor))
    for i in 1:length(species)
        plot!(p, [0.0], [0.0], lab=specieslabels[i],
              c=speciescolor[species[i]],
              leg=:topright, legendfontsize=5,
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


