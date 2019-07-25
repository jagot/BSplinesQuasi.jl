using ArnoldiMethod
using LinearAlgebra

function test_bspline_derivatives(t, f::Function, g::Function, h::Function)
    R = BSpline(t)
    D = Derivative(axes(R,1))

    ∇ = R'D*R
    ∇² = R'D'D*R

    S = R'R
    S⁻¹ = qr(S)

    fv = R \ f
    gv = similar(fv)
    hv = similar(fv)

    mul!(gv, ∇, fv)
    ldiv!(S⁻¹, gv)
    mul!(hv, ∇², fv)
    ldiv!(S⁻¹, hv)

    δg = gv - (R \ g)
    δh = hv - (R \ h)

    fv,gv,hv,δg,δh
end

function error_slope(loghs,ϵ)
    # To avoid the effect of round-off errors on the order
    # estimation.
    i = argmin(abs.(ϵ)) - 1

    ([loghs[1:i] ones(i)] \ log10.(abs.(ϵ[1:i])))[1]
end

function compute_derivative_errors(a, b, k, Ns, f::Function, g::Function, h::Function)
    errors = map(Ns) do N
        t = LinearKnotSet(k, a, b, N)
        fv,gv,hv,δg,δh = test_bspline_derivatives(t, f, g, h)

        [norm(δg) norm(δh)]
    end |> e -> vcat(e...)

    ϵg = errors[:,1]
    ϵh = errors[:,2]

    loghs = log10.(1.0 ./ Ns)
    pg = error_slope(loghs, ϵg)
    ph = error_slope(loghs, ϵh)

    ϵg,ϵh,pg,ph
end

function derivative_test_functions(d)
    a,b = √d*[-1,1]

    # The functions need to vanish at the boundaries, for the
    # derivative approximation to be valid (only Dirichlet0 boundary
    # conditions implemented).
    f = x -> exp(-1/(d-x^2))
    g = x -> -2*exp(-1/(d-x^2))*x/((d-x^2)^2)
    h = x -> -2*exp(-1/(d-x^2))*(d^2 + 2*(d-1)*x^2-3x^4)/((d-x^2)^4)

    f,g,h,a,b
end
