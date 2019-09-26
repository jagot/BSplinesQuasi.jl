# * Diagonal operators

function Matrix(f::Function, L::BSpline{T}, R::BSpline{T}) where T
    χ = L.B
    ξ = Diagonal(f.(R.x))*R.B
    k = max(order(L),order(R))
    m,n = size(L,2),size(R,2)

    S = BandedMatrix(Zeros{T}(m,n), (k-1,k-1))
    overlap_matrix!(S, χ, ξ, weights(R))
end

function Matrix(f::Function,
                L::RestrictedQuasiArray{<:Any,2,<:BSpline},
                R::RestrictedQuasiArray{<:Any,2,<:BSpline})
    L′,Lrestriction = L.args
    R′,Rrestriction = R.args
    La,Lb = restriction_extents(Lrestriction)
    Ra,Rb = restriction_extents(Rrestriction)
    M = Matrix(f, L′, R′)
    M[1+La:end-Lb,1+Ra:end-Rb]
end

"""
    Matrix(f, B)

Generate the matrix representing the diagonal operator `f` in the
B-spline space `B`.

# Examples

```jldoctest
julia> Matrix(x -> x, BSpline(LinearKnotSet(3, 0, 1, 3)))
5×5 BandedMatrices.BandedMatrix{Float64,Array{Float64,2},Base.OneTo{Int64}}:
 0.0037037    0.00462963  0.000925926   ⋅           ⋅
 0.00462963   0.0259259   0.0236111    0.00138889   ⋅
 0.000925926  0.0236111   0.0916667    0.0458333   0.00462963
  ⋅           0.00138889  0.0458333    0.0851852   0.0342593
  ⋅            ⋅          0.00462963   0.0342593   0.062963
```
"""
function Matrix(f::Function, B::BSpline{T}) where T
    χ = B.B
    ξ = Diagonal(f.(B.x))*B.B
    k = order(B)
    n = size(B,2)

    S = BandedMatrix(Zeros{T}(n,n), (k-1,k-1))
    overlap_matrix!(S, χ, ξ, weights(B))
end

function Matrix(f::Function, B::RestrictedQuasiArray{<:Any,2,<:BSpline})
    B′,restriction = B.args
    a,b = restriction_extents(restriction)
    M = Matrix(f, B′)
    M[1+a:end-b,1+a:end-b]
end
