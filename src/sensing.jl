export
    sense_f!, sense_∇f!, sense_f_∇f!

function sense_f!(bacterium::B, f::F, affect!; kwargs...) where {B<:AbstractBacterium,F<:AbstractAnalyticField}
    ϕ = f.ϕ(bacterium; kwargs...)
    affect!(bacterium, ϕ)
end # function

function sense_∇f!(bacterium::B, f::F, affect!; kwargs...) where {B<:AbstractBacterium,F<:AbstractAnalyticField}
    ∇ϕ = f.∇ϕ(bacterium; kwargs...)
    affect!(bacterium, ∇ϕ)
end # function

function sense_f_∇f!(bacterium::B, f::F, affect!; kwargs...) where {B<:AbstractBacterium,F<:AbstractAnalyticField}
    ϕ = f.ϕ(bacterium; kwargs...)
    ∇ϕ = f.∇ϕ(bacterium; kwargs...)
    affect!(bacterium, ϕ, ∇ϕ)
end # function
