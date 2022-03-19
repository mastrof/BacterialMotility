export
    step!

@doc raw"""
    step!(bacterium::B, dt, f::F) where {B<:AbstractBacterium,F<:AbstractField}
    step!(bacteria::T, dt, f::F) where {B<:AbstractBacterium,T<:AbstractVector{B},F<:AbstractField}

Perform a single integration step over the interval dt, composed as:
    1) run
    2) sensing
    3) reorientation attempt
"""
function step!(bacterium::B, dt, f::F) where {B<:AbstractBacterium,F<:AbstractField}
    bacterium.run!(bacterium, dt)
    bacterium.sense!(bacterium, f)
    if rand() < bacterium.state["ReorientationRate"] * dt
        bacterium.turn!(bacterium)
    end # if
end # function


function step!(bacteria::T, dt, f::F) where {B<:AbstractBacterium,T<:AbstractVector{B},F<:AbstractField}
    for bacterium in bacteria
        step!(bacterium, dt, f)
    end # for
end # function

