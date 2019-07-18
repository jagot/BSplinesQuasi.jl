const RestrictionMatrix = BandedMatrix{<:Int, <:FillArrays.Ones}

const RestrictedBasis{B<:Basis} = Mul{<:Any,<:Tuple{B, <:RestrictionMatrix}}
const AdjointRestrictedBasis{B<:Basis} = Mul{<:Any,<:Tuple{<:Adjoint{<:Any,<:RestrictionMatrix}, <:QuasiAdjoint{<:Any,B}}}

const RestrictedQuasiArray{T,N,B<:Basis} = MulQuasiArray{T,N,<:RestrictedBasis{B}}
const AdjointRestrictedQuasiArray{T,N,B<:Basis} = MulQuasiArray{T,N,<:AdjointRestrictedBasis{B}}

const BasisOrRestricted{B<:Basis} = Union{B,RestrictedBasis{<:B},<:RestrictedQuasiArray{<:Any,<:Any,<:B}}
const AdjointBasisOrRestricted{B<:Basis} = Union{<:QuasiAdjoint{<:Any,B},AdjointRestrictedBasis{<:B},<:AdjointRestrictedQuasiArray{<:Any,<:Any,<:B}}

unrestricted_basis(R::AbstractQuasiMatrix) = R
unrestricted_basis(R::RestrictedBasis) = first(R.args)
unrestricted_basis(R::RestrictedQuasiArray) = unrestricted_basis(R.applied)

function restriction_extents(restriction::RestrictionMatrix)
    a = restriction.l
    b = restriction.raxis.stop - restriction.data.axes[2].stop - a
    a,b
end
