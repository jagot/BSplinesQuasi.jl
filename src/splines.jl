"""
    BSpline(t, x, w, B, S)

Basis structure for the B-splines generated by the knot set `t`. `x`
and `w` are the associated quadrature roots and weights, respectively,
the columns of `B` correspond to the B-splines resolved on the
quadrature roots, and `S` is the (banded) B-spline overlap matrix.
"""
struct BSpline{T,R<:Real,
               KnotSet<:AbstractKnotSet{<:Any,<:Any,<:Any,R},
               XV<:AbstractVector{T},
               WV<:AbstractVector{R},
               BM<:AbstractMatrix{T},
               SM<:AbstractMatrix{T}} <: Basis{T}
    t::KnotSet
    x::XV
    w::WV
    B::BM
    S::SM
end

function BSpline(t::AbstractKnotSet{k}, x::AbstractVector{T}, w::AbstractVector) where {k,T}
    nf = numfunctions(t)
    B = spzeros(T, length(x), nf)

    nei = nonempty_intervals(t)
    # This is the amount of quadrature points per interval (assuming
    # the same amount per interval, which is currently the case), this
    # could conceivably be passed as an argument, instead.
    N = findlast(e -> e < t[nei[1]+1], x)

    ix = 1
    for i ∈ nei
        # These are the functions which are non-zero on the interval
        # t[i]..t[i+1].
        for j ∈ max(1,i-k+1):min(i,nf)
            eⱼ = UnitVector{T}(nf, j)
            for l ∈ ix:ix+N-1
                B[l,j] = deBoor(t, eⱼ, x[l], i)
            end
        end
        ix += N
    end
    # For basis function i, loop over intervals j ∈ i:i+k
    # Find all quadrature roots within interval j
    # Evaluate B_{i,k}(x_j)
    S = BandedMatrix(Zeros{T}(nf, nf), (k-1,k-1))
    for i ∈ 1:nf
        for j ∈ max(i-k+1,1):min(i+k-1,nf)
            S[i,j] = B[:,i]' * Diagonal(w) * B[:,j]
        end
    end
    BSpline(t, x, w, B, S)
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
    N = numfunctions(B′.t)
    show(io, B′)
    write(io, ", restricted to basis functions $(1+a)..$(N-b) $(a>0 || b>0 ? "⊂" : "⊆") 1..$(N)")
end

restriction_extents(B::BSpline) = 0,0
restriction_extents(B::RestrictedQuasiArray{<:Any,2,<:BSpline}) =
    restriction_extents(B.applied.args[2])

locs(B::BSpline) = B.x

locs(B::RestrictedQuasiArray{<:Any,2,<:BSpline}) =
    first(B.applied.args).x

weights(B::BSpline) = B.w

weights(B::RestrictedQuasiArray{<:Any,2,<:BSpline}) =
    first(B.applied.args).w

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

@inline function Base.getindex(B::RestrictedBasis{<:BSpline{T}}, x::Real, sel::AbstractVector) where {T}
    B′,restriction = B.args
    B′[x,sel .+ restriction.l]
end

function getindex(B::BSpline{T}, x::AbstractRange, j::Integer) where T
    χ = spzeros(T, length(x))
    basis_function!(χ, B, x, j)
    χ
end

@inline function Base.getindex(B::RestrictedBasis{<:BSpline{T}}, x::Real, j::Integer) where {T}
    B′,restriction = B.args
    B′[x,j+restriction.l]
end

getindex(B::BSplineOrRestricted, x, ::Colon) =
    getindex(B, x, axes(B,2))

# * Types

const SplineArray{T,N,B<:BSplineOrRestricted} = MulQuasiArray{T,N,<:Mul{<:Any,<:Tuple{B,<:AbstractArray{T,N}}}}
const SplineVector{T,B<:BSplineOrRestricted} = SplineArray{T,1,B}
const SplineMatrix{T,B<:BSplineOrRestricted} = SplineArray{T,2,B}
const SplineVecOrMat{T,B<:BSplineOrRestricted} = Union{SplineVector{T,B},SplineMatrix{T,B}}

Base.show(io::IO, spline::SplineVector) =
    write(io, "Spline on $(spline.applied.args[1])")

Base.show(io::IO, spline::SplineMatrix) =
    write(io, "$(size(spline, 2))d spline on $(spline.applied.args[1])")

# * Function interpolation

Base.:(\ )(B::BSpline, f::Function) = B.B \ f.(B.x)
function Base.:(\ )(B::RestrictedBSpline, f::Function)
    # B′,restriction = B.applied.args
    # a,b = restriction_extents(restriction)
    x = locs(B)
    V = B[x,:]
    V \ f.(x)
end

export BSpline
