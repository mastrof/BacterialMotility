export sense!

#=
function sense!(bacterium::B, f::F, affect!) where {B<:AbstractBacterium,F<:AbstractAnalyticField}
    ϕ = f.ϕ(bacterium)
    ∇ϕ = f.∇ϕ(bacterium)
    ∂ₜϕ = f.∂ₜϕ(bacterium)
    affect!(bacterium, ϕ, ∇ϕ, ∂ₜϕ)
end # function
=#

function sense!(bs, i, affect!)
    ϕ = bs.f.ϕ(bs, i)
    ∇ϕ = bs.f.∇ϕ(bs, i)
    ∂ₜϕ = bs.f.∂ₜϕ(bs, i)
    affect!(bs.population[i], ϕ, ∇ϕ, ∂ₜϕ)
end # function
