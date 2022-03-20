export
    sense_f!, sense_∇f!, sense_f_∇f!

function sense_f!(bacterium::B, f::F, affect!; kwargs...) where {B<:AbstractBacterium,F<:AbstractAnalyticField}
    ϕ = f.ϕ(bacterium)
    affect!(bacterium, ϕ; kwargs...)
end # function

function sense_∇f!(bacterium::B, f::F, affect!; kwargs...) where {B<:AbstractBacterium,F<:AbstractAnalyticField}
    ∇ϕ = f.∇ϕ(bacterium)
    affect!(bacterium, ∇ϕ; kwargs...)
end # function

function sense_f_∇f!(bacterium::B, f::F, affect!; kwargs...) where {B<:AbstractBacterium,F<:AbstractAnalyticField}
    ϕ = f.ϕ(bacterium)
    ∇ϕ = f.∇ϕ(bacterium)
    affect!(bacterium, ϕ, ∇ϕ; kwargs...)
end # function
