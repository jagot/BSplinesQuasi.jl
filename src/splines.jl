struct BSpline{T,R<:Real,
               KnotSet<:AbstractKnotSet{<:Any,<:Any,<:Any,R},
               XV<:AbstractVector{T},
               WV<:AbstractVector{R}} <: Basis{T}
    t::KnotSet
    x::XV
    w::WV
end

"""
    BSpline(t[, k′=3])

Create the B-spline basis corresponding to the knot set `t`. `k′` is
the highest polynomial order of operators for which it should be
possible to compute the matrix elements exactly (via Gauß–Legendre
quadrature). The default `k′=3` corresponds to operators O(x²).
"""
BSpline(t::AbstractKnotSet, k′::Integer=3) = BSpline(t, lgwt(t, k′)...)

const RestrictedBSpline{T} = Union{RestrictedBasis{<:BSpline{T}},<:RestrictedQuasiArray{<:Any,<:Any,<:BSpline{T}}}
const AdjointRestrictedBSpline{T} = Union{AdjointRestrictedBasis{<:BSpline{T}},<:AdjointRestrictedQuasiArray{<:Any,<:Any,<:BSpline{T}}}

const BSplineOrRestricted{T} = BasisOrRestricted{<:BSpline{T}}
const AdjointBSplineOrRestricted{T} = AdjointBasisOrRestricted{<:BSpline{T}}

# * Properties

axes(B::BSpline) = (Inclusion(first(B.t)..last(B.t)), Base.OneTo(numfunctions(B.t)))
size(B::BSpline) = (ℵ₁, numfunctions(B.t))
size(B::RestrictedQuasiArray{<:Any,2,<:BSpline}) = (ℵ₁, length(B.applied.args[2].data))
==(A::BSpline,B::BSpline) = A.t == B.t
==(A::BSplineOrRestricted,B::BSplineOrRestricted) = unrestricted_basis(A) == unrestricted_basis(B)

order(B::BSpline) = order(B.t)

function show(io::IO, B::BSpline{T}) where T
    write(io, "BSpline{$(T)} basis with $(B.t)")
end

function show(io::IO, B::RestrictedQuasiArray{T,2,BSpline{T}}) where T
    B′,restriction = B.applied.args
    a,b = restriction_extents(restriction)
    N = length(B′.x)
    show(io, B′)
    write(io, ", restricted to basis functions $(1+a)..$(N-b) $(a>0 || b>0 ? "⊂" : "⊆") 1..$(N)")
end

restriction_extents(B::BSpline) = 0,0
restriction_extents(B::RestrictedQuasiArray{<:Any,2,<:BSpline}) =
    restriction_extents(B.applied.args[2])

locs(B::BSpline) = B.x

function locs(B::RestrictedQuasiArray{<:Any,2,<:BSpline})
    B′,restriction = B.applied.args
    a,b = BSplineQuasi.restriction_extents(restriction)
    B′.x[1+a:end-b]
end

IntervalSets.leftendpoint(B::BSpline) = B.x[1]
IntervalSets.rightendpoint(B::BSpline) = B.x[end]

IntervalSets.leftendpoint(B::RestrictedQuasiArray{<:Any,2,<:BSpline}) =
    leftendpoint(B.applied.args[1])
IntervalSets.rightendpoint(B::RestrictedQuasiArray{<:Any,2,<:BSpline}) =
    rightendpoint(B.applied.args[1])

# # * Basis functions

"""
    deBoor(t, c, x[, k])

Evaluate the spline given by the knot set `t` and the set of control
points `c` at `x` using de Boor's algorithm. `k` is the index of the
knot interval containing `x`.
"""
function deBoor(t::AbstractKnotSet, c::AbstractVector, x,
                i=find_interval(t, x))
    isnothing(i) && return zero(eltype(t))
    k = order(t)
    nc = length(c)
    k == 1 && return nc < i ? 0 : c[i]

    α = [r > 0 && r ≤ nc ? c[r] : zero(eltype(c))
         for r ∈ i-k+1:i]
    nt = length(t)
    nf = numfunctions(t)
    for j = 1:k-1
        for r = i:-1:max(i-k+j,1)
            jjj = r+k-j
            (jjj > nt || jjj < 1) && continue
            r′ = r - i + k
            r ≠ 1 && r′ == 1 && continue

            a = t[r+k-j]
            b = t[r]
            
            α[r′] = if r == 1
                (b-x)*α[r′]/(b-a)
            elseif r == nf + j
                (x-a)*α[r′-1]/(b-a)
            else
                ((x-a)*α[r′-1] + (b-x)*α[r′])/(b-a)
            end
        end
    end

    α[end]
end

getindex(B::BSpline{T}, x::Real, j::Integer) where T =
    deBoor(B.t, UnitVector{T}(size(B,2), j),
           x, find_interval(B.t, x))

function basis_function!(χ, B::BSpline{T}, x::AbstractRange, j) where T
    eⱼ = UnitVector{T}(size(B,2), j)
    for (is,k) ∈ within_support(x, B.t, j)
        for i in is
            χ[i] = deBoor(B.t, eⱼ, x[i], k)
        end
    end
end

function getindex(B::BSpline{T}, x::AbstractRange, sel::AbstractVector) where T
    χ = spzeros(T, length(x), length(sel))
    for j in sel
        basis_function!(view(χ, :, j), B, x, j)
    end
    χ
end

function getindex(B::BSpline{T}, x::AbstractRange, j::Integer) where T
    χ = spzeros(T, length(x))
    basis_function!(χ, B, x, j)
    χ
end

getindex(B::BSpline{T}, x, ::Colon) where T =
    getindex(B, x, axes(B,2))

export BSpline
