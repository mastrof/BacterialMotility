export
    AbstractField, EmptyField,
    AbstractAnalyticField, AnalyticField,
    Concentration_SteadyDiffusionSphericalSource_3D

abstract type AbstractField end

abstract type AbstractAnalyticField <: AbstractField end

@with_kw struct AnalyticField <: AbstractAnalyticField
    ϕ::Function = zeroscalar
    ∇ϕ::Function = zerovector
    ∂ₜϕ::Function = zeroscalar
end # struct

EmptyField() = AnalyticField()



@inline function ϕ_steadydiffspherical(bs, i; kwargs...)
    b = bs.population[i]
    C = kwargs[:C]
    Cbg = kwargs[:Cbg]
    R = kwargs[:R]
    r = norm(b.r)
    Cbg + (C-Cbg) * R / r
end # function
                  
@inline function ∇ϕ_steadydiffspherical(bs, i; kwargs...)
    b = bs.population[i]
    C = kwargs[:C]
    Cbg = kwargs[:Cbg]
    R = kwargs[:R]
    r = norm(b.r)
    r³ = r*r*r
    -C*R / r³ .* b.r
end # function

@with_kw struct Concentration_SteadyDiffusionSphericalSource_3D <: AbstractAnalyticField
    R::Float64 # μm, source radius
    C::Float64 # μM, concentration at source surface
    Cbg::Float64 = 0.0 # μM, background concentration
    ϕ::Function = (bs,i) -> ϕ_steadydiffspherical(bs, i; C=C, Cbg=Cbg, R=R) # concentration field
    ∇ϕ::Function = (bs,i) -> ∇ϕ_steadydiffspherical(bs, i; C=C, Cbg=Cbg, R=R)  # concentration gradient (radial towards center)
    ∂ₜϕ::Function = zeroscalar # no time derivative
end # struct


#=
abstract type AbstractNumericField <: AbstractField

@with_kw struct NumericField <: AbstractNumericField
    ϕ
    ∇ϕ
    ∂ₜϕ
end # struct
=#
