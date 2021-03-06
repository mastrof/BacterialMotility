export
    BacterialSystem,
    integrate!,
    step!

@doc raw"""
    struct BacterialSystem <: AbstractBacterialSystem

Collects all the information about the system to be integrated.
    - `clock`: tracks the number of integration steps
    - `stop`: can be turned to true with callbacks to stop integration when a certain condition is satisfied
    - `field`: external field(s) in which bacteria move (e.g. a concentration gradient) (<:AbstractField)
    - `boundary_conditions!`: function that enforces boundary conditions on the bacteria; takes only single bacterium
    - `callback_inner`: user-defined function called after the timestep of each bacterium; takes BacterialSystem and index of current bacterium (Int)
    - `callback_outer`: user-defined function called after the whole population has undergone one timestep; takes only BacterialSystem
    - `user_parameters`: collection of mutables for user-defined routines (e.g. Dict, user-defined struct, ...)
    - `population`: AbstractVector of <:AbstractBacterium whose motion is to be integrated
"""
@with_kw struct BacterialSystem <: AbstractBacterialSystem
    clock::Vector{Int} = [0]
    stop::Vector{Bool} = [false]
    field = EmptyField() # AbstractField
    boundary_conditions! = dummy_boundary_conditions # Function
    callback_inner = dummy_callback_inner # Function
    callback_outer = dummy_callback_outer # Function
    user_parameters = Dict() # might use better name
    population # AbstractVector{AbstractBacterium}
end # struct


@doc raw"""
    step!(bs::AbstractBacterialSystem, i::Int)

Perform a single integration step of bacterium indicated by index `i` in bacterial system `bs`.
"""
function step!(bs::AbstractBacterialSystem, i::Int)
    bacterium = bs.population[i]
    dt = bacterium.state["IntegrationTimestep"]
    bacterium.run!(bacterium, dt)
    bs.boundary_conditions!(bacterium)
    bacterium.sense!(bs, i)
    ν = bacterium.state["ReorientationRate"]
    if rand() < ν * dt
        bacterium.turn!(bacterium)
    end # if
    bs.callback_inner(bs, i)
    rotational_diffusion!(bacterium)
end # function


@doc raw"""
    integrate!(bs::AbstractBacterialSystem, nsteps::Int)

Integrate the bacterial system `bs` for `nsteps` steps.
"""
function integrate!(bs::AbstractBacterialSystem, nsteps::Int)
    for n in 1:nsteps
        for i in eachindex(bs.population)
            step!(bs, i)
        end # for
        bs.clock[1] += 1
        bs.callback_outer(bs)
        bs.stop[1] && break # if satisfied, interrupt integration
    end # for
end # function


#== NOW BROKEN, NEEDS REVIEW ==#
#=
function integrate!(bs::BacterialSystem, T::Float64)
    num_bacteria = length(bs.population)
    Δt = [b.state["IntegrationTimestep"] for b in bs.population]
    nsteps = T ./ Δt # steps for each bacterium to reach time T
    N = round(Int, maximum(nsteps))
    δ = nsteps ./ N
    accumulator = zeros(num_bacteria)
    for n in 1:N
        for i in 1:num_bacteria
            accumulator[i] += δ[i]
            if accumulator[i] ≥ 1
                step!(bs.population[i], bs.field;
                      bcs! = bs.boundary_conditions!, callback=bs.callback_inner, kwargs...)
                accumulator[i] -= 1
            end # if
        end # for
        bs.clock[1] += 1
        bs.callback_outer(bs; kwargs...)
    end # for
end # function
=#
