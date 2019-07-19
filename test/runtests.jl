using BSplinesQuasi
using Test
using LinearAlgebra
using SparseArrays
using Compat

function vecdist(a::AbstractVector, b::AbstractVector,
                 ϵ = eps(eltype(a)))
    δ = √(sum(abs2, a-b))
    δ, δ/√(sum(abs2, a .+ ϵ))
end

include("knot_sets.jl")
include("splines.jl")
