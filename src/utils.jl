#=
  Utilities used throughout the library, not exported
=#

dummy(args...; kwargs...) = nothing
dummy_run(b::AbstractBacterium, dt) = nothing
dummy_turn(b::AbstractBacterium) = nothing
dummy_boundary_conditions(b::AbstractBacterium) = nothing
dummy_sense(bs::AbstractBacterialSystem, i) = nothing
dummy_callback_inner(bs::AbstractBacterialSystem, i::Int) = nothing
dummy_callback_outer(bs::AbstractBacterialSystem) = nothing

zerofunc(args...; kwargs...) = 0.0
zeroscalar(bs::AbstractBacterialSystem, i::Int) = 0.0
zerovector(bs::AbstractBacterialSystem, i::Int) = zeros(size(first(bs.population).r))

struct Degenerate{T<:Real} <: ContinuousUnivariateDistribution
    x::T
end # struct

Base.rand(d::Degenerate) = d.x
