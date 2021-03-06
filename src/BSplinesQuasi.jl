module BSplinesQuasi

function __init__()
    @warn "The BSplinesQuasi.jl package has been deprecated in favour of JuliaApproximation/CompactBases.jl"
    nothing
end

import Base: axes, size, ==, getindex, checkbounds, copyto!, similar, diff, show
import Base.Broadcast: materialize

using ContinuumArrays
import ContinuumArrays: Basis, ℵ₁, Derivative, @simplify

using QuasiArrays
import QuasiArrays: AbstractQuasiMatrix, QuasiAdjoint, MulQuasiArray, Inclusion, ApplyQuasiArray

using BandedMatrices
using SparseArrays

using IntervalSets

using LazyArrays
import LazyArrays: ⋆
using FillArrays

using LinearAlgebra
import LinearAlgebra: Matrix, dot

using FastGaussQuadrature

using RecipesBase
using Compat

"""
    UnitVector{T}(N, k)

Helper vector type of length `N` where the `k`th element is `one(T)`
and all the others `zero(T)`.
"""
struct UnitVector{T} <: AbstractVector{T}
    N::Int
    k::Int
end

Base.size(e::UnitVector) = (e.N,)
Base.getindex(e::UnitVector{T}, i::Int) where T = i == e.k ? one(T) : zero(T)
Base.show(io::IO, e::UnitVector{T}) where T = write(io, "ê{$T}($(e.k)|$(e.N))")

include("knot_sets.jl")
include("quadrature.jl")
include("restricted_bases.jl")
include("splines.jl")
include("operators.jl")
include("derivatives.jl")

export Derivative

end # module
