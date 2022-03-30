export
    step!, multistep!

@doc raw"""
    step!(bacterium::B, f::F; callback) where {B<:AbstractBacterium,F<:AbstractField}
    step!(population::T, f::F; callback_inner, callback_outer) where {B<:AbstractBacterium,T<:AbstractVector{B},F<:AbstractField}

Perform a single integration step, composed as:
    1) run
    2) sensing
    3) reorientation attempt
    4) callback
"""
function step!(bacterium::B, f::F=EmptyField(); callback=dummy) where {B<:AbstractBacterium,F<:AbstractField}
    dt = bacterium.state["IntegrationTimestep"]
    bacterium.run!(bacterium, dt)
    bacterium.sense!(bacterium, f)
    ν = bacterium.state["ReorientationRate"]
    if rand() < ν * dt
        bacterium.turn!(bacterium)
    end # if
    callback(bacterium, f)
end # function

function step!(population::T, f::F=EmptyField(); callback_inner=dummy, callback_outer=dummy) where {B<:AbstractBacterium,T<:AbstractVector{B},F<:AbstractField}
    for bacterium in population
        step!(bacterium, f; callback=callback_inner)
    end # for
    callback_outer(population, f)
end # function


@doc raw"""
    multistep!(population::V, T, f::F=EmptyField; callback_inner, callback_outer) where {B<:AbstractBacterium,V<:AbstractVector{B},F<:AbstractField}

Integrate the equations of motion of the `population` in field `f`, up to time `T` with multi-timestep integration.
It is recommended (but not required) to have IntegrationTimesteps that are divisors of `T` for maximum accuracy.
The integration time `T` should be given in the same units of the IntegrationTimesteps.
If IntegrationTimesteps are all equal, the `step!` function should be used instead (although this will not lead to errors, only lower performance).
"""
function multistep!(population::V, T, f::F=EmptyField(); callback_inner=dummy, callback_outer=dummy) where {B<:AbstractBacterium,V<:AbstractVector{B},F<:AbstractField}
    Nb = length(population)
    Δt = [b.state["IntegrationTimestep"] for b in population]
    nsteps = T ./ Δt # steps required for each bacterium to reach time T
    N = maximum(nsteps)
    δ = nsteps ./ N
    accumulator = zeros(Nb)
    for n in 1:N
        for i in 1:Nb
            accumulator[i] += δ[i]
            if accumulator[i] ≥ 1
                step!(population[i], f; callback=callback_inner)
                accumulator[i] -= 1
            end # if
        end # for
        callback_outer()
    end # for
end # function
