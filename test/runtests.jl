using BSplinesQuasi
using Test
using LinearAlgebra
using BandedMatrices
using SparseArrays
using Compat

import ContinuumArrays: ℵ₁

function vecdist(a::AbstractVector, b::AbstractVector,
                 ϵ = eps(eltype(a)))
    δ = √(sum(abs2, a-b))
    δ, δ/√(sum(abs2, a .+ ϵ))
end

include("knot_sets.jl")
include("splines.jl")
include("derivatives.jl")
