#=
  Utilities used throughout the library, not exported
=#

dummy(args...; kwargs...) = nothing

zerofunc(args...; kwargs...) = 0.0

struct Degenerate{T<:Real} <: ContinuousUnivariateDistribution
    x::T
end # struct

Base.rand(d::Degenerate) = d.x
