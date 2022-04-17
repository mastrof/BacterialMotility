export sense!

function sense!(bacterium::B, f::F, affect!; kwargs...) where {B<:AbstractBacterium,F<:AbstractAnalyticField}
    ϕ = f.ϕ(bacterium; kwargs...)
    ∇ϕ = f.∇ϕ(bacterium; kwargs...)
    ∂ₜϕ = f.∂ₜϕ(bacterium; kwargs...)
    affect!(bacterium, ϕ, ∇ϕ, ∂ₜϕ)
end # function
