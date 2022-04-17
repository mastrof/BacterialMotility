export
    AbstractField, EmptyField,
    AbstractAnalyticField, AnalyticField,
    Concentration_SteadyDiffusionSphericalSource_3D

abstract type AbstractField end

abstract type AbstractAnalyticField <: AbstractField end

@with_kw struct AnalyticField <: AbstractAnalyticField
    ϕ::Function = zerofunc
    ∇ϕ::Function = zerofunc
    ∂ₜϕ::Function = zerofunc
end # struct

EmptyField() = AnalyticField()

@with_kw struct Concentration_SteadyDiffusionSphericalSource_3D <: AbstractAnalyticField
    R::Float64 # μm, source radius
    C::Float64 # μM, concentration at source surface
    Cbg::Float64 = 0.0 # μM, background concentration
    ϕ::Function = (b; kwargs...) -> Cbg + (C-Cbg)*R/norm(b.r) # concentration field
    ∇ϕ::Function = (b; kwargs...) -> -C*R/(norm(b.r)^2) # concentration gradient (radial towards center)
    ∂ₜϕ::Function = zerofunc # no time derivative
end # struct


#=
abstract type AbstractNumericField <: AbstractField

@with_kw struct NumericField <: AbstractNumericField
    ϕ
    ∇ϕ
    ∂ₜϕ
end # struct
=#
