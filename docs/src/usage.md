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

## Splines

A spline is constructed as a linear combination of B-splines:

$$\begin{equation}
s(x) = \sum_{j=1}^{n_{tk}} \B{j}{k}(x)c_j \defd B\vec{c}.
\end{equation}$$

This is easily done as
```julia
julia> c = sin.(1:size(B,2))
8-element Array{Float64,1}:
  0.8414709848078965
  0.9092974268256817
  0.1411200080598672
 -0.7568024953079282
 -0.9589242746631385
 -0.27941549819892586
  0.6569865987187891
  0.9893582466233818

julia> s = B*c
ContinuumArrays.QuasiArrays.ApplyQuasiArray{Float64,1,LazyArrays.Applied{ContinuumArrays.BasisStyle,typeof(*),Tuple{BSpline{Float64,Float64,LinearKnotSet{4,4,4,Float64,StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},Array{Float64,1},Array{Float64,1}},Array{Float64,1}}}}(BSpline{Float64} basis with LinearKnotSet(Float64) of order k = 4 (cubic) on 0.0..1.0 (5 intervals)⋆[0.841471, 0.909297, 0.14112, -0.756802, -0.958924, -0.279415, 0.656987, 0.989358])
```

Naturally, we can evaluate the spline similarly to above:
```julia
julia> s[0.3]
-0.28804656969083225

julia> s[0.3:0.1:0.8]
6-element Array{Float64,1}:
 -0.28804656969083225
 -0.6408357079724976
 -0.8250002333223664
 -0.8119858486932346
 -0.5856964504991202
 -0.15856643671353235
```

If many different spline sharing the same B-splines are going to be
evaluated, it is usually more efficient to evaluate the basis
functions once and reuse them:

```julia
julia> χ = B[0:0.25:1.0, :]
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

julia> χ*c
5-element Array{Float64,1}:
  0.8414709848078965
 -0.06366510061255656
 -0.8250002333223664
 -0.39601358159504896
  0.9893582466233818
```

![One-dimensional spline](figures/spline-1d.svg)

It is then trivial to extend this to two dimensions:

```julia
julia> c′ = [sin.(1:size(B,2)) tan.(1:size(B,2))]
8×2 Array{Float64,2}:
  0.841471   0.540302
  0.909297  -0.416147
  0.14112   -0.989992
 -0.756802  -0.653644
 -0.958924   0.283662
 -0.279415   0.96017
  0.656987   0.753902
  0.989358  -0.1455

julia> s′ = B*c′
ContinuumArrays.QuasiArrays.ApplyQuasiArray{Float64,2,LazyArrays.Applied{ContinuumArrays.BasisStyle,typeof(*),Tuple{BSpline{Float64,Float64,LinearKnotSet{4,4,4,Float64,StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},Array{Float64,1},Array{Float64,1}},Array{Float64,2}}}}(BSpline{Float64} basis with LinearKnotSet(Float64) of order k = 4 (cubic) on 0.0..1.0 (5 intervals)⋆[0.841471 0.540302; 0.909297 -0.416147; … ; 0.656987 0.753902; 0.989358 -0.1455])

julia> size(s′)
(ContinuumArrays.AlephInfinity{1}(), 2)
```

![Two-dimensional spline](figures/spline-2d.svg)