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


#== custom motile patterns ==#
custom_tumble!(bacterium) = tumble!(bacterium, Uniform(-π/2, π/2))

#== population setup ==#
Δt = 0.1 # s
D_rot = 0.1 # rad²/s
s1 = properties("IntegrationTimestep" => Δt,
                "RotationalDiffusivity" => D_rot)
s2 = properties("IntegrationTimestep" => Δt,
                "RotationalDiffusivity" => D_rot)
s3 = properties("IntegrationTimestep" => Δt,
                "RotationalDiffusivity" => D_rot,
                "ReorientationRate" => 0.5)

species = ["RT_1", "RT_2", "RRF_1"]
N = 20
pop1 = [Bacterium{d}(
    id = species[1], r = randposition(), v = randvelocity(36.0),
    run! = run!, turn! = tumble!, state = copy(s1)) for _ in 1:N]
pop2 = [Bacterium{d}(
    id = species[2], r = randposition(), v = randvelocity(25.0),
    run! = run!, turn! = custom_tumble!, state = copy(s2)) for _ in 1:2N]
pop3 = [Bacterium{d}(
    id = species[3], r = randposition(), v = randvelocity(45.0),
    run! = run!, turn! = reverse_flick!, state = copy(s3)) for _ in 1:N]

population = vcat(pop1, pop2, pop3)
num_bacteria = length(population)




#== callbacks ==#
function pbc!(b, Lmin, Lmax)
    ΔL = Lmax - Lmin
    for i in 1:length(b.r)
        if b.r[i] < Lmin
            b.r[i] += ΔL
        elseif b.r[i] > Lmax
            b.r[i] -= ΔL
        end # if
    end # for
end # function

pbc!(b) = pbc!(b, -L/2, L/2)

function callback_inner(b,f; kwargs...)
    rotational_diffusion!(b)
    #pbc!(b)
end # function

function hunt!(bs, radius_perception, prey, hunter)
    bacteria_prey = [b for b in bs.population if b.id == prey]
    bacteria_hunter = [b for b in bs.population if b.id == hunter]
    nb_prey = length(bacteria_prey)
    nb_hunter = length(bacteria_hunter)
    for i in 1:nb_hunter
        r_i = bacteria_hunter[i].r
        relative_distances = [norm(r_i .- bacteria_prey[j].r) for j in 1:nb_prey]
        Δr, j = findmin(relative_distances)
        Δr > radius_perception && continue
        U_i = norm(bacteria_hunter[i].v)
        r_j = bacteria_prey[j].r
        unitvector_ij = (r_j .- r_i) ./ Δr
        perception_noise = rand(Normal(0,1/3), length(unitvector_ij))
        unitvector_ij .= (unitvector_ij .+ perception_noise) ./ norm(unitvector_ij .+ perception_noise)
        bacteria_hunter[i].v .= unitvector_ij .* U_i
        rotational_diffusion!(bacteria_hunter[i]) # reapply rotational diffusion which gets overwritten by interaction
        #println(Δr, " ", i, " ", j, " ", unitvector_ij, " ", U_i)
    end # for
end # function

function save_traj!(traj, bs)
    nsteps, num_bacteria, d = size(traj)
    t = bs.clock[1]
    for i in 1:num_bacteria
        traj[t,i,:] .= bs.population[i].r
    end # for
end # function

RPs = [30.0, 50.0] # μm
nsteps = 600
trajectories = zeros(nsteps, num_bacteria, d)

function callback_outer(bs, RPs, preys, hunters; kwargs...)
    for (RP,p,h) in zip(RPs, preys, hunters)
        hunt!(bs, RP, p, h)
    end # for
    save_traj!(trajectories, bs)
end # function
callback_outer(bs; kwargs...) = callback_outer(bs, RPs, ["RT_1", "RRF_1"], ["RT_2", "RT_1"]; kwargs...)


#== simulation ==#
bs = BacterialSystem(callback_inner = callback_inner,
                     callback_outer = callback_outer,
                     population = population)

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
#specieslabels = ["Fast", "Slow", "Slower"]
linecolor = [speciescolor[bact.id] for _ in 1:1, bact in population]

bgcolor = RGB(0.07, 0.07, 0.07)

#== plot animation ==#
ltail = 40
for t in 2:2:nsteps
    p = plot(;plot_style(:Dark2)...)
    plot!(p, lims=(-L/2,L/2), aspect_ratio=1, axis=false,
          bgcolor=bgcolor, size=(600,600))
    t0 = max(t-ltail, 1)
    tailplot!(p, trajectories, t0, t; lc=linecolor)
    headplot!(p, trajectories, t; mc=permutedims(linecolor))
    #=
    for i in 1:length(species)
        plot!(p, [0.0], [0.0], lab=specieslabels[i],
              lc=speciescolor[species[i]],
              legend=:topright, legendfontsize=5,
              legendfontcolor=:white)
    end # for
    =#
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
