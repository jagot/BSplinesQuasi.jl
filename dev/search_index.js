var documenterSearchIndex = {"docs":
[{"location":"splines/#Spline-creation-and-evaluation-1","page":"Spline creation & evaluation","title":"Spline creation & evaluation","text":"","category":"section"},{"location":"splines/#","page":"Spline creation & evaluation","title":"Spline creation & evaluation","text":"A spline is constructed as a linear combination of B-splines:","category":"page"},{"location":"splines/#","page":"Spline creation & evaluation","title":"Spline creation & evaluation","text":"beginequation\nlabeleqnspline\ns(x) = sum_j=1^n_tk Bjk(x)c_j defd Bvecc\nendequation","category":"page"},{"location":"splines/#","page":"Spline creation & evaluation","title":"Spline creation & evaluation","text":"This is easily done as","category":"page"},{"location":"splines/#","page":"Spline creation & evaluation","title":"Spline creation & evaluation","text":"julia> k = 4\n4\n\njulia> t = LinearKnotSet(k, 0, 1, 5)\n12-element LinearKnotSet{4,4,4,Float64,StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}}:\n 0.0\n 0.0\n 0.0\n 0.0\n 0.2\n 0.4\n 0.6\n 0.8\n 1.0\n 1.0\n 1.0\n 1.0\n\njulia> B = BSpline(t)\nBSpline{Float64} basis with LinearKnotSet(Float64) of order k = 4 (cubic) on 0.0..1.0 (5 intervals)\n\njulia> size(B)\n(ContinuumArrays.AlephInfinity{1}(), 8)\n\njulia> c = sin.(1:size(B,2))\n8-element Array{Float64,1}:\n  0.8414709848078965\n  0.9092974268256817\n  0.1411200080598672\n -0.7568024953079282\n -0.9589242746631385\n -0.27941549819892586\n  0.6569865987187891\n  0.9893582466233818\n\njulia> s = B*c\nSpline on BSpline{Float64} basis with LinearKnotSet(Float64) of order k = 4 (cubic) on 0.0..1.0 (5 intervals)","category":"page"},{"location":"splines/#","page":"Spline creation & evaluation","title":"Spline creation & evaluation","text":"Naturally, we can evaluate the spline similarly to above:","category":"page"},{"location":"splines/#","page":"Spline creation & evaluation","title":"Spline creation & evaluation","text":"julia> s[0.3]\n-0.28804656969083225\n\njulia> s[0.3:0.1:0.8]\n6-element Array{Float64,1}:\n -0.28804656969083225\n -0.6408357079724976\n -0.8250002333223664\n -0.8119858486932346\n -0.5856964504991202\n -0.15856643671353235","category":"page"},{"location":"splines/#","page":"Spline creation & evaluation","title":"Spline creation & evaluation","text":"If many different splines sharing the same B-splines are going to be evaluated, it is usually more efficient to evaluate the basis functions once and reuse them:","category":"page"},{"location":"splines/#","page":"Spline creation & evaluation","title":"Spline creation & evaluation","text":"julia> χ = B[0:0.25:1.0, :]\n5×8 SparseArrays.SparseMatrixCSC{Float64,Int64} with 14 stored entries:\n  [1, 1]  =  1.0\n  [2, 2]  =  0.105469\n  [2, 3]  =  0.576823\n  [3, 3]  =  0.0208333\n  [2, 4]  =  0.315104\n  [3, 4]  =  0.479167\n  [4, 4]  =  0.00260417\n  [2, 5]  =  0.00260417\n  [3, 5]  =  0.479167\n  [4, 5]  =  0.315104\n  [3, 6]  =  0.0208333\n  [4, 6]  =  0.576823\n  [4, 7]  =  0.105469\n  [5, 8]  =  1.0\n\njulia> χ*c\n5-element Array{Float64,1}:\n  0.8414709848078965\n -0.06366510061255656\n -0.8250002333223664\n -0.39601358159504896\n  0.9893582466233818","category":"page"},{"location":"splines/#","page":"Spline creation & evaluation","title":"Spline creation & evaluation","text":"(Image: One-dimensional spline)","category":"page"},{"location":"splines/#","page":"Spline creation & evaluation","title":"Spline creation & evaluation","text":"It is then trivial to extend this to two dimensions:","category":"page"},{"location":"splines/#","page":"Spline creation & evaluation","title":"Spline creation & evaluation","text":"julia> c̃ = [sin.(1:size(B,2)) tan.(1:size(B,2))]\n8×2 Array{Float64,2}:\n  0.841471   1.55741\n  0.909297  -2.18504\n  0.14112   -0.142547\n -0.756802   1.15782\n -0.958924  -3.38052\n -0.279415  -0.291006\n  0.656987   0.871448\n  0.989358  -6.79971\n\njulia> s̃ = B*c̃\n2d spline on BSpline{Float64} basis with LinearKnotSet(Float64) of order k = 4 (cubic) on 0.0..1.0 (5 intervals)\n\njulia> size(s̃)\n(ContinuumArrays.AlephInfinity{1}(), 2)","category":"page"},{"location":"splines/#","page":"Spline creation & evaluation","title":"Spline creation & evaluation","text":"(Image: Two-dimensional spline)","category":"page"},{"location":"theory/#Theory-1","page":"Theory","title":"Theory","text":"","category":"section"},{"location":"theory/#","page":"Theory","title":"Theory","text":"The underlying equations for the definition and efficient evaluation of B-splines are introduced. Equations labelled (dB.A.ij) refer to equation (ij) in chapter A of","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"Carl de Boor (2001). A Practical Guide to Splines. New York: Springer. ISBN: 978-0-387-95366-3.","category":"page"},{"location":"theory/#Definitions-1","page":"Theory","title":"Definitions","text":"","category":"section"},{"location":"theory/#","page":"Theory","title":"Theory","text":"The order of the polynomial is designated k, equal to the degree + 1, i.e. a parabolic polynomial is of order 3.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"The B-splines can be defined through the Cox–de Boor recursion relation: beginequation tagdBIX14 Bjk defd omega_jkBjk-1 + (1-omega_j+1k)Bj+1k-1 endequation where beginequation tagdBIX15 omega_jk(x) defd fracx-t_jt_j+k-1-t_j endequation and beginequation tagdBIX11 Br1(x) = begincases 1  xin t_rt_r+1)\n0  textrmelse\nendcasesquad rin1n_t-1 endequation","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"Given a knot vector vect of length n_t, there are n_tkdefd n_t-k functions of order k. This implies that there is a highest order k a given knot set vector can support, i.e. k_textrmmax = n_t - 1.","category":"page"},{"location":"theory/#Knot-sets-1","page":"Theory","title":"Knot sets","text":"","category":"section"},{"location":"theory/#","page":"Theory","title":"Theory","text":"As seen above, the B-splines are completely defined by the knot set vector vect.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"The knot set vect=bmat123456 gives rise to the following cardinal splines, of orders k=15: (Image: Cardinal splines)","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"By increasing the multiplicity of some knots, the continuity of the splines can be controlled. E.g. the knot set vect=bmat01134666, will yield the following splines: (Image: Discontinuous splines)","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"Lastly, it is very common to pad the knot set such that the first and last knot have multiplicity k; this simplifies the implementation of boundary conditions when B-splines are used for solving differential equations: (Image: Full multiplicity splines)","category":"page"},{"location":"theory/#de-Boor's-algorithm-1","page":"Theory","title":"de Boor's algorithm","text":"","category":"section"},{"location":"theory/#","page":"Theory","title":"Theory","text":"An efficient way of evaluating the B-splines is given by de Boor's algorithm. The algorithm described in the Wikipedia article [1] assumes full multiplicity at the endpoints of the knot set, i.e. that the first and last points are repeated k times. In the description of the same algorithm, de Boor additionally assumes a strictly increasing knot set, i.e. t_i+1t_i, do avoid any divisions by zero. Since one of the goals for this package is maximum flexibility in choosing the knot set, a derivation of an only slightly more general version of de Boor's algorithm follows.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"Given a knot set vect of length n_t, the spline f(x) is given by beginequation tagdBX23⅓  ⅔ beginaligned f(x) = sum_r=1^n_t-k alpha_r Brk(x)\n= sum_r=1^n_t-k alpha_rfracx-t_rt_r+k-1-t_rBrk-1(x)+ sum_r=1^n_t-k alpha_rfract_r+k-xt_r+k-t_r+1Br+1k-1(x)\n= sum_r=1^n_t-k alpha_rfracx-t_rt_r+k-1-t_rBrk-1(x)+ sum_r=2^n_t-k+1 alpha_r-1fract_r+k-1-xt_r+k-1-t_rBrk-1(x)\n= fracx-t_1t_k-t_1alpha_1B1k-1(x) + left sum_r=2^n_t-k frac(x-t_r)alpha_r+(t_r+k-1-x)alpha_r-1t_r+k-1-t_rBrk-1(x) right + fract_n_t-xt_n_t-t_n_t-k+1alpha_n_t-kBn_t-k+1k-1(x)\n= sum_i=1^n_t-k+1 alpha_i^2(x) Bik-1(x) endaligned endequation where beginequation tagdBX24 alpha_r^2(x) defd begincases displaystyle fracx-t_1t_k-t_1alpha_1  r=12ex displaystyle frac(x-t_r)alpha_r+(t_r+k-1-x)alpha_r-1t_r+k-1-t_r  rin2n_t-k2ex displaystyle fract_n_t-xt_n_t-t_n_t-k+1alpha_n_t-k  r = n_t - k + 1 endcases endequation We have thus reexpressed the spline function f(x) of order k as a linear combination of B-splines of order k-1. We can generalize this, to reexpress f(x) as a linear combination of B-splines of order k-j, with expansion coefficients beginequation tagdBX26 alpha_r^j+1(x) defd begincases displaystyle fracx-t_1t_1+k-j-t_1alpha_1^j(x)  r=12ex displaystyle frac(x-t_r)alpha_r^j(x)+(t_r+k-j-x)alpha_r-1^j(x)t_r+k-j-t_r  rin2n_t-k+j-12ex displaystyle fract_n_t-xt_n_t-t_n_t-k+jalpha_n_t-k+j-1^j(x)  r = n_t - k + j endcases endequation","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"The difference between this derivation and those in the Wikipedia article and de Boor (2001), is that we here explicitly consider the limits of the sum imposed by the length of the knot set and the order k; this introduces the special cases for r=1n_t-k+j.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"An important philosophical difference between the Cox–de Boor recursion relation and de Boor's algorithm, is that whereas the former is a linear combination of basis functions evaluated at certain position x, the latter is linear combination of intervals (since the first-order functions Bj1 are non-zero within one interval only, and they are mutual orthogonal), with polynomial expansion coefficients alpha_i^k(x). To evaluate the spline function f(x), we first find the interval i which contains x. Even if the knot set is only non-decreasing, i.e. not strictly increasing, the interval containing x is uniquely defined, since there is only one for which t_i leq x  t_i+1; if the knot t_i has a multiplicity higher than unity, the additional intervals cannot contain x, since they are empty: t_i-1leq x  t_i = varnothing if t_i-1 = t_i. By finding the last i, for which t_i leq x, we thus guarantee that no divisions by zero will occur.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"[1]: NB that the Wikipedia article uses 0-based indexing, whereas de Boor and BSplinesQuasi.jl use 1-based indexing.","category":"page"},{"location":"theory/#Integrals-1","page":"Theory","title":"Integrals","text":"","category":"section"},{"location":"theory/#","page":"Theory","title":"Theory","text":"Since the B-splines are piecewise polynomials, they can be exactly integrated using Gauß–Legendre quadrature; an N-point quadrature can integrate a polynomial of degree 2N-1 exactly. We are usually interested in integrals on the form","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"beginequation\nmatrixelBj_1koperatorABj_2kequiv\nintdiffx\nconjBj_1k(x)\noperatorA\nBj_2k(x)\nendequation","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"where operatorA is some operator. Assuming that operatorA can be approximated by a polynomial of order k, we need an N-point quadrature, such that 2N-1geq 2(k-1)+(k-1), e.g. for operatorAsim x^2, we choose N=ceilfrac2k+12=k+1.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"For every non-empty interval generated by the knot set, we setup a Gauß–Legendre quadrature, such that an integral is approximated as","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"beginequation\nintlimits_a^b diffxf(x)approx\nfracb-a2sum_i=1^n w_i\nfleft(fracb-a2x_i+fraca+b2right)\nendequation","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"within each interval.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"If we again consider the knot set vect=bmat01134666 and allow operators of maximal polynomial order k=3, we get the following distribution quadrature points:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"(Image: Quadrature points)","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"Note that no quadrature points were generated for the intervals t_it_i+1), i=267, since those intervals are empty.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"With the quadrature in place, it becomes very easy to compute the overlap matrix:","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"beginequation\nmatS_ijk defd braketBikBjk\napprox sum_l w_l conjBik(x_l) Bjk(x_l)\nendequation","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"For the knot set above, we find","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"matS =\nbmat06        0222222    00444444   00         00\n 0222222   0466667    0307407    00037037   00\n 00444444   0307407    0962963    0307407    00444444\n 00        00037037   0307407    0466667    0222222\n 00        00         00444444   0222222    04","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"from which we see that the individual B-splines Bik are non-zero only on the interval t_it_i+k), except for the last B-spline that is non-zero also at the end of the interval, t_it_i+k.","category":"page"},{"location":"theory/#","page":"Theory","title":"Theory","text":"If we want to employ two B-spline sets of different orders, we must make sure they share the same knot set and quadrature points (and that the latter support the combined polynomial order).","category":"page"},{"location":"usage/#Usage-1","page":"Basis creation","title":"Usage","text":"","category":"section"},{"location":"usage/#","page":"Basis creation","title":"Basis creation","text":"We first load the package","category":"page"},{"location":"usage/#","page":"Basis creation","title":"Basis creation","text":"using BSplinesQuasi","category":"page"},{"location":"usage/#","page":"Basis creation","title":"Basis creation","text":"Then we define a linearly spaced knot set between 0 and 1 with five intervals and cubic splines. By default, full multiplicity of the endpoints is assumed.","category":"page"},{"location":"usage/#","page":"Basis creation","title":"Basis creation","text":"julia> t = LinearKnotSet(4, 0, 1, 5)\n12-element LinearKnotSet{4,4,4,Float64,StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}}:\n 0.0\n 0.0\n 0.0\n 0.0\n 0.2\n 0.4\n 0.6\n 0.8\n 1.0\n 1.0\n 1.0\n 1.0\n\njulia> B = BSpline(t)\nBSpline{Float64} basis with LinearKnotSet(Float64) of order k = 4 (cubic) on 0.0..1.0 (5 intervals)\n\njulia> size(B)\n(ContinuumArrays.AlephInfinity{1}(), 8)","category":"page"},{"location":"usage/#","page":"Basis creation","title":"Basis creation","text":"The last statement means that B is a quasimatrix with a continuous first dimension which is spanned by 8 basis functions.","category":"page"},{"location":"usage/#Evaluation-of-B-splines-1","page":"Basis creation","title":"Evaluation of B-splines","text":"","category":"section"},{"location":"usage/#","page":"Basis creation","title":"Basis creation","text":"It is the possible to query the values of e.g. the first basis functions at some values of x:","category":"page"},{"location":"usage/#","page":"Basis creation","title":"Basis creation","text":"julia> B[0:0.1:1,1]\n11-element SparseArrays.SparseVector{Float64,Int64} with 2 stored entries:\n  [1 ]  =  1.0\n  [2 ]  =  0.125","category":"page"},{"location":"usage/#","page":"Basis creation","title":"Basis creation","text":"We can also evaluate all the B-splines at the same time:","category":"page"},{"location":"usage/#","page":"Basis creation","title":"Basis creation","text":"julia> B[0:0.25:1,:]\n5×8 SparseArrays.SparseMatrixCSC{Float64,Int64} with 14 stored entries:\n  [1, 1]  =  1.0\n  [2, 2]  =  0.105469\n  [2, 3]  =  0.576823\n  [3, 3]  =  0.0208333\n  [2, 4]  =  0.315104\n  [3, 4]  =  0.479167\n  [4, 4]  =  0.00260417\n  [2, 5]  =  0.00260417\n  [3, 5]  =  0.479167\n  [4, 5]  =  0.315104\n  [3, 6]  =  0.0208333\n  [4, 6]  =  0.576823\n  [4, 7]  =  0.105469\n  [5, 8]  =  1.0","category":"page"},{"location":"usage/#","page":"Basis creation","title":"Basis creation","text":"Since the B-splines have compact support, they are only locally non-zero, hence the sparse storage.","category":"page"},{"location":"usage/#","page":"Basis creation","title":"Basis creation","text":"Finally, we can compute the value of a single B-spline, a range, or all of them at a single point x:","category":"page"},{"location":"usage/#","page":"Basis creation","title":"Basis creation","text":"julia> B[0.5,4]\n0.47916666666666674\n\njulia> B[0.5,:]\n8-element Array{Float64,1}:\n 0.0\n 0.0\n 0.020833333333333325\n 0.47916666666666674\n 0.4791666666666666\n 0.020833333333333322\n 0.0\n 0.0\n\njulia> B[0.5,4:8]\n5-element Array{Float64,1}:\n 0.47916666666666674\n 0.4791666666666666\n 0.020833333333333322\n 0.0\n 0.0","category":"page"},{"location":"function_approximation/#Function-approximation-1","page":"Function approximation","title":"Function approximation","text":"","category":"section"},{"location":"function_approximation/#","page":"Function approximation","title":"Function approximation","text":"To approximate a mathematical function on a B-spline basis, we can simply solve for the coefficients:","category":"page"},{"location":"function_approximation/#","page":"Function approximation","title":"Function approximation","text":"julia> t = LinearKnotSet(7, 0, 7, 10);\n\njulia> B = BSpline(t)\nBSpline{Float64} basis with LinearKnotSet(Float64) of order k = 7 on 0.0..7.0 (10 intervals)\n\njulia> c = B \\ sin\n16-element Array{Float64,1}:\n  3.674947319163227e-7\n  0.1166649763895669\n  0.3500043164970867\n  0.6828451145424481\n  1.0237038168861594\n  1.1358262863822566\n  0.7361804983706153\n -0.009705277902046651\n -0.7510239693627963\n -1.139127153806961\n -0.991477377385712\n -0.4798514317407342\n  0.024153163433647172\n  0.3716550472969573\n  0.5690327347258564\n  0.6569863188695337","category":"page"},{"location":"function_approximation/#","page":"Function approximation","title":"Function approximation","text":"The top panel shows the expansion coefficents and the reconstructed function, the middle panel the reconstruction error, and the bottom panel the underlying basis functions.","category":"page"},{"location":"function_approximation/#","page":"Function approximation","title":"Function approximation","text":"(Image: Function interpolation by B-splines)","category":"page"},{"location":"function_approximation/#","page":"Function approximation","title":"Function approximation","text":"Since the sine function is non-zero at x=7, it is important that our basis set includes a B-spline that supports this, hence the full multiplicity of last knot. The sine function is zero at the first knot, however, something that is reflected in the fact that the first expansion coefficient is almost zero. In problems where vanishing boundary conditions are stipulated, this can be enforced by dropping the first/last spline:","category":"page"},{"location":"function_approximation/#","page":"Function approximation","title":"Function approximation","text":"julia> t = LinearKnotSet(7, 0.0, 1.0, 6);\n\njulia> B = BSpline(t,3)\nBSpline{Float64} basis with LinearKnotSet(Float64) of order k = 7 on 0.0..1.0 (6 intervals)\n\njulia> B̃ = B[:,2:end-1]\nBSpline{Float64} basis with LinearKnotSet(Float64) of order k = 7 on 0.0..1.0 (6 intervals), restricted to basis functions 2..11 ⊂ 1..12","category":"page"},{"location":"function_approximation/#","page":"Function approximation","title":"Function approximation","text":"We can now compare how well the restricted basis can reconstruct different functions, compared to the unrestricted one:","category":"page"},{"location":"function_approximation/#","page":"Function approximation","title":"Function approximation","text":"julia> f1 = x -> sin(2π*x)\n#121 (generic function with 1 method)\n\njulia> f2 = x -> cos(2π*x)\n#123 (generic function with 1 method)\n\njulia> c1 = B \\ f1\n12-element Array{Float64,1}:\n  8.595306123332097e-6\n  0.17449541085301523\n  0.5236843390613484\n  0.9897088684171409\n  1.265948043821677\n  0.6905122508504199\n -0.6905122508504197\n -1.2659480438216766\n -0.9897088684171416\n -0.5236843390613481\n -0.1744954108530153\n -8.59530612374404e-6\n\njulia> c2 = B \\ f2\n12-element Array{Float64,1}:\n  0.9999984763007934\n  1.000010902252537\n  0.9268416666707764\n  0.5980612634315096\n -0.1993704588377225\n -1.1959917207779527\n -1.195991720777952\n -0.1993704588377244\n  0.5980612634315112\n  0.9268416666707747\n  1.000010902252537\n  0.9999984763007937\n\njulia> c̃1 = B̃ \\ f1\n10-element Array{Float64,1}:\n  0.17450648562854001\n  0.5236722006849034\n  0.9897209715549453\n  1.2659362792307434\n  0.6905240028400118\n -0.690524002840012\n -1.265936279230742\n -0.989720971554948\n -0.5236722006849018\n -0.1745064856285418\n\njulia> c̃2 = B̃ \\ f2\n10-element Array{Float64,1}:\n  2.1869627421644986\n -0.18872692489766923\n  1.4520927974960849\n -0.6709944899983972\n -1.057228429296141\n -1.057228429296146\n -0.6709944899983946\n  1.4520927974960847\n -0.1887269248976696\n  2.1869627421644977","category":"page"},{"location":"function_approximation/#","page":"Function approximation","title":"Function approximation","text":"(Image: Reconstruction of function interpolated on restricted bases)","category":"page"},{"location":"function_approximation/#","page":"Function approximation","title":"Function approximation","text":"As is to be expected, the sine function is perfectly reconstructed in both cases, whilst the cosine fails spectactularly in the restricted case.","category":"page"},{"location":"#BSplinesQuasi.jl-1","page":"Home","title":"BSplinesQuasi.jl","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"BSplinesQuasi is a package for B-splines in the framework of ContinuumArrays","category":"page"},{"location":"#","page":"Home","title":"Home","text":"(Image: BSplinesQuasi logo)","category":"page"},{"location":"#Index-of-types-and-functions-1","page":"Home","title":"Index of types and functions","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Modules = [BSplinesQuasi]","category":"page"}]
}
