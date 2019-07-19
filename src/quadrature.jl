using FastGaussQuadrature

#=
Gaußian quadrature:

\[\int\limits_a^b dx\;f(x)\approx
\frac{b-a}{2}\sum_{i=1}^n w_i
f\left(\frac{b-a}{2}x_i+\frac{a+b}{2}\right)\]
=#
function lgwt!(x::AbstractVector{T}, w::AbstractVector{T},
               xs::AbstractVector{T}, ws::AbstractVector{T},
               a::T=zero(T), b::T=one(T)) where {T<:Real}
    xs .= 0.5((b-a)*x .+ a .+ b)
    ws .= 0.5(b-a)*w
    nothing
end

"""
    num_quadrature_points(k, k′)

The number of quadrature points to compute the matrix element of an
operator of polynomial order `k′` with respect to a basis of order
`k`.
"""
function num_quadrature_points(k, k′)
    N2 = 2*(k-1) + k′
    N2>>1 + N2&1
end

function lgwt(t::AbstractKnotSet{k,ml,mr,T}, k′) where {k,ml,mr,T}
    N = num_quadrature_points(k, k′)
    x, w = gausslegendre(N)

    nei = nonempty_intervals(t)
    ni = length(nei)
    xo = zeros(T, ni*length(x))
    wo = zeros(T, ni*length(x))

    for (i,j) in enumerate(nei)
        sel = (i-1)*N+1 : i*N
        lgwt!(x, w,
              view(xo, sel), view(wo, sel),
              t[j], t[j+1])
    end

    xo,wo
end
