# * Derivative types

# ** First derivatives

const FlatFirstDerivative = Mul{<:Any, <:Tuple{
    <:QuasiAdjoint{<:Any, <:BSpline},
    <:Derivative,
    <:BSpline}}
const FlatRestrictedFirstDerivative = Mul{<:Any, <:Tuple{
    <:Adjoint{<:Any,<:RestrictionMatrix},
    <:QuasiAdjoint{<:Any, <:BSpline},
    <:Derivative,
    <:BSpline,
    <:RestrictionMatrix}}

const LazyFirstDerivative = Mul{<:Any, <:Tuple{
    <:Mul{<:Any, <:Tuple{
        <:QuasiAdjoint{<:Any, <:BSpline},
        <:Derivative}},
    <:BSpline}}

const LazyRestrictedFirstDerivative = Mul{<:Any, <:Tuple{
    <:Mul{<:Any,<:Tuple{
        <:MulQuasiArray{<:Any, 2, <:Mul{<:Any, <:Tuple{
            <:Adjoint{<:Any,<:RestrictionMatrix},
            <:QuasiAdjoint{<:Any,<:BSpline}}}},
        <:Derivative}},
    <:RestrictedQuasiArray{<:Any,2,<:BSpline}}}

const FirstDerivative = Union{FlatFirstDerivative, FlatRestrictedFirstDerivative,
                              LazyFirstDerivative, LazyRestrictedFirstDerivative}

# ** Second derivatives

const FlatSecondDerivative = Mul{<:Any, <:Tuple{
    <:QuasiAdjoint{<:Any, <:BSpline},
    <:QuasiAdjoint{<:Any, <:Derivative},
    <:Derivative,
    <:BSpline}}
const FlatRestrictedSecondDerivative = Mul{<:Any, <:Tuple{
    <:Adjoint{<:Any,<:RestrictionMatrix},
    <:QuasiAdjoint{<:Any, <:BSpline},
    <:QuasiAdjoint{<:Any, <:Derivative},
    <:Derivative,
    <:BSpline,
    <:RestrictionMatrix}}

const LazySecondDerivative = Mul{<:Any, <:Tuple{
    <:Mul{<:Any, <:Tuple{
        <:Mul{<:Any, <:Tuple{
            <:QuasiAdjoint{<:Any, <:BSpline}, <:QuasiAdjoint{<:Any, <:Derivative}}},
        <:Derivative}},
    <:BSpline}}

const LazyRestrictedSecondDerivative = Mul{<:Any, <:Tuple{
    <:Mul{<:Any,<:Tuple{
        <:Mul{<:Any,<:Tuple{
            <:MulQuasiArray{<:Any, 2, <:Mul{<:Any, <:Tuple{
                <:Adjoint{<:Any,<:RestrictionMatrix},
                <:QuasiAdjoint{<:Any,<:BSpline}}}},
            <:QuasiAdjoint{<:Any,<:Derivative}}},
        <:Derivative}},
    <:RestrictedQuasiArray{<:Any,2,<:BSpline}}}

const SecondDerivative = Union{FlatSecondDerivative,FlatRestrictedSecondDerivative,
                               LazySecondDerivative,LazyRestrictedSecondDerivative}

const FirstOrSecondDerivative = Union{FirstDerivative,SecondDerivative}

difforder(::FirstDerivative) = 1
difforder(::SecondDerivative) = 2

# * Materialization

function diffop!(dest::BandedMatrix, B::BSplineOrRestricted, o)
    k = order(B)
    (bandwidth(dest,1) ≤ k-1) && (bandwidth(dest,2) ≤ k-1) ||
        throw(DimensionMismatch("Insufficient bandwidths of destination matrix"))

    d,r = divrem(o,2)
    a,b = d,d+r

    # ∂' = -∂ ⇒ for weak Laplacians, we negate the left basis.
    χ = (iseven(a) ? 1 : -1)*basis_functions(B, a)
    ξ = basis_functions(B, b)

    overlap_matrix!(dest, χ, ξ, weights(B))

    dest
end

function copyto!(dest::BandedMatrix, M::FirstOrSecondDerivative)
    axes(dest) == axes(M) || throw(DimensionMismatch("axes must be same"))
    B = basis(M)

    diffop!(dest, basis(M), difforder(M))
end

basis(M::Union{FlatFirstDerivative,LazyFirstDerivative,
               FlatSecondDerivative,LazySecondDerivative}) = last(M.args)

basis(M::Union{FlatRestrictedFirstDerivative, FlatRestrictedSecondDerivative}) =
    M.args[end-1]*M.args[end]

function similar(M::FirstOrSecondDerivative, ::Type{T}) where T
    B = basis(M)
    k = order(B)
    BandedMatrix(Zeros{T}(size(M)), (k-1,k-1))
end
materialize(M::FirstOrSecondDerivative) = copyto!(similar(M, eltype(M)), M)
