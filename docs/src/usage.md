# Usage

We first load the package
```julia
using BSplinesQuasi
```

Then we define a linearly spaced knot set between `0` and `1` with
five intervals and cubic splines. By default, full multiplicity of the
endpoints is assumed.
```julia
julia> t = LinearKnotSet(4, 0, 1, 5)
12-element LinearKnotSet{4,4,4,Float64,StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}}:
 0.0
 0.0
 0.0
 0.0
 0.2
 0.4
 0.6
 0.8
 1.0
 1.0
 1.0
 1.0

julia> B = BSpline(t)
BSpline{Float64} basis with LinearKnotSet(Float64) of order k = 4 (cubic) on 0.0..1.0 (5 intervals)

julia> size(B)
(ContinuumArrays.AlephInfinity{1}(), 8)
```
The last statement means that `B` is a quasimatrix with a continuous
first dimension which is spanned by 8 basis functions.

## Evaluation of B-splines

It is the possible to query the values of e.g. the first basis functions at some values of `x`:

```julia
julia> B[0:0.1:1,1]
11-element SparseArrays.SparseVector{Float64,Int64} with 2 stored entries:
  [1 ]  =  1.0
  [2 ]  =  0.125
```

We can also evaluate all the B-splines at the same time:

```julia
julia> B[0:0.25:1,:]
5×8 SparseArrays.SparseMatrixCSC{Float64,Int64} with 14 stored entries:
  [1, 1]  =  1.0
  [2, 2]  =  0.105469
  [2, 3]  =  0.576823
  [3, 3]  =  0.0208333
  [2, 4]  =  0.315104
  [3, 4]  =  0.479167
  [4, 4]  =  0.00260417
  [2, 5]  =  0.00260417
  [3, 5]  =  0.479167
  [4, 5]  =  0.315104
  [3, 6]  =  0.0208333
  [4, 6]  =  0.576823
  [4, 7]  =  0.105469
  [5, 8]  =  1.0
```

Since the B-splines have compact support, they are only locally
non-zero, hence the sparse storage.

Finally, we can compute the value of a single B-spline, a range, or
all of them at a single point `x`:

```julia
julia> B[0.5,4]
0.47916666666666674

julia> B[0.5,:]
8-element Array{Float64,1}:
 0.0
 0.0
 0.020833333333333325
 0.47916666666666674
 0.4791666666666666
 0.020833333333333322
 0.0
 0.0

julia> B[0.5,4:8]
5-element Array{Float64,1}:
 0.47916666666666674
 0.4791666666666666
 0.020833333333333322
 0.0
 0.0
```