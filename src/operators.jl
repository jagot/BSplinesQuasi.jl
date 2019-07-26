# * Diagonal operators

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
    B′,restriction = B.applied.args
    a,b = restriction_extents(restriction)
    M = Matrix(f, B′)
    M[1+a:end-b,1+a:end-b]
end
