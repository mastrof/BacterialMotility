#=
  Utilities used throughout the library, not exported
=#

dummy(args...; kwargs...) = nothing
dummy_run(b::AbstractBacterium, dt) = nothing
dummy_turn(b::AbstractBacterium) = nothing
dummy_sense(bs::BacterialSystem, i) = nothing
dummy_callback_inner(bs::BacterialSystem, i::Int) = nothing
dummy_callback_outer(bs::BacterialSystem) = nothing
dummy_boundary_conditions(b::AbstractBacterium) = nothing

zerofunc(args...; kwargs...) = 0.0
zeroscalar(bs::BacterialSystem, i::Int) = 0.0
zerovector(bs::BacterialSystem, i::Int) = zeros(size(first(bs.population).r))

struct Degenerate{T<:Real} <: ContinuousUnivariateDistribution
    x::T
end # struct

Base.rand(d::Degenerate) = d.x
