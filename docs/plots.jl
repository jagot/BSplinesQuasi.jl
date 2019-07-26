using PyPlot
using Jagot.plotting
plot_style("ggplot")
using PyPlotRecipes
using PyCall
PathEffects = pyimport("matplotlib.patheffects")
using Statistics
using Random
using Colors

using LinearAlgebra
using BSplinesQuasi

using ArnoldiMethod

function mean_color(color::String)
    c = parse(Colorant, color)
    mean([c.r,c.g,c.b])
end

lerp(a,b,t) = (1-t)*a + t*b

mean_position(x, ϕ) = ϕ'*Diagonal(x)*ϕ/(ϕ'ϕ)

function cfigure(fun::Function, figname; clear=true, tight=true, kwargs...)
    figure(figname; kwargs...)
    clear && clf()
    fun()
    tight && tight_layout()
end

function csubplot(fun::Function, args...; nox=false, noy=false)
    ax = subplot(args...)
    fun()
    if nox
        no_tick_labels(:x)
        xlabel("")
    end
    if noy
        no_tick_labels(:y)
        ylabel("")
    end
    ax
end


function Bspline_text(x, ϕ, j, k, color)
    xⱼ = mean_position(x, ϕ)
    txt = text(xⱼ,0.7maximum(ϕ),
               latexstring("\\mathrm{B}_{$j,$k}"),
               horizontalalignment="center",
               color=color)
    gray = 1.0 - mean_color(color)
    outline = RGB(lerp(gray, (gray<0.5 ? 0 : 1), 0.6)*[1,1,1]...)
    txt.set_path_effects([PathEffects.Stroke(linewidth=1,
                                             foreground="#$(hex(outline))"),
                          PathEffects.Normal()])
end

function cardinal_splines()
    a,b = 1.0,6.0
    x = range(a, stop=b, length=1001)
    figure("cardinal splines",figsize=(7,9))
    clf()
    t = 0
    for k = 1:5
        t = ArbitraryKnotSet(k, a:b, 1, 1)
        R = BSpline(t)
        χ = R[x,:]
        subplot(5,1,k)
        for j = 1:size(χ,2)
            ϕ = view(χ, :, j)
            l=plot(x, ϕ)[1]
            Bspline_text(x, ϕ, j, k, l.get_color())
        end
        ylabel(latexstring("k = $k"))
        margins(0.1,0.1)
        no_tick_labels(:x)
    end
    xlabel(L"x")
    tight_layout()
    savefig("docs/src/figures/cardinal-splines.svg")
end

function discontinuous_splines()
    t = ArbitraryKnotSet(3, [0.0, 1, 1, 3, 4, 6], 1, 3)
    r = range(first(t), stop=last(t), length=301)
    R = BSpline(t)
    χ = R[r, :]

    figure("basis functions",figsize=(7,6))
    clf()
    subplot(211)
    plot(r, χ, length(r) < 30 ? ".-" : "-")
    margins(0.1, 0.1)
    no_tick_labels()
    subplot(212)
    rplot(t)
    margins(0.1, 0.1)
    xlabel(L"r")
    ylabel("Multiplicity")
    tight_layout()
    savefig("docs/src/figures/discontinuous-splines.svg")
end

function full_multiplicity_splines()
    a,b = 1.0,6.0
    x = range(a, stop=b, length=1001)
    figure("splines",figsize=(7,9))
    clf()
    t = 0
    for k = 1:5
        t = ArbitraryKnotSet(k, a:b)
        R = BSpline(t)
        χ = R[x,:]
        subplot(5,1,k)
        for j = 1:size(χ,2)
            ϕ = view(χ, :, j)
            l=plot(x, ϕ)[1]
            Bspline_text(x, ϕ, j, k, l.get_color())
        end
        ylabel(latexstring("k = $k"))
        margins(0.1,0.1)
        no_tick_labels(:x)
    end
    xlabel(L"x")
    tight_layout()
    savefig("docs/src/figures/full-multiplicity-splines.svg")
end

function spline1d()
    k = 4
    t = LinearKnotSet(k, 0, 1, 5)
    B = BSpline(t)
    x = range(0, stop=1, length=301)
    χ = B[x, :]

    i = 1:size(B,2)
    c = sin.(i)
    s = χ*c

    figure("spline 1d", figsize=(7,9))
    clf()
    subplot(311)
    plot(x, s)
    plot([mean_position(x, view(χ, :, j)) for j in i], c, "s-")
    no_tick_labels()
    ylabel(L"s(x)")
    margins(0.1,0.1)
    subplot(312)
    for j = 1:size(χ,2)
        ϕ = view(χ, :, j)
        l=plot(x, ϕ)[1]
        Bspline_text(x, ϕ, j, k, l.get_color())
    end
    margins(0.1,0.1)
    no_tick_labels()
    subplot(313)
    rplot(t)
    margins(0.1,0.1)
    xlabel(L"x")
    tight_layout()

    savefig("docs/src/figures/spline-1d.svg")
end

function spline2d()
    k = 4
    t = LinearKnotSet(k, 0, 1, 5)
    B = BSpline(t)
    x = range(0, stop=1, length=301)
    χ = B[x, :]

    i = 1:size(B,2)
    c = [sin.(i) tan.(i)]
    s = χ*c

    figure("spline 2d", figsize=(6,6))
    clf()
    plot(s[:,1], s[:,2])
    plot(c[:,1], c[:,2], "s-")
    xlabel(L"x")
    ylabel(L"y")
    margins(0.1,0.1)
    tight_layout()

    savefig("docs/src/figures/spline-2d.svg")
end

function logo()
    t = ArbitraryKnotSet(3, [0.0, 1, 1, 3, 4, 6], 1, 3)
    r = range(first(t), stop=last(t), length=301)
    R = BSpline(t)
    χ = R[r, :]

    figure("logo",figsize=(7,3))
    clf()
    plot(r, χ, length(r) < 30 ? ".-" : "-", linewidth=4)
    margins(0.1, 0.1)
    axis("off")
    tight_layout()
    savefig("docs/src/assets/logo.svg", transparent=true)
end

function quadrature_points()
    t = ArbitraryKnotSet(3, [0.0, 1, 1, 3, 4, 6], 1, 3)
    x = range(first(t), stop=last(t), length=301)
    B = BSpline(t,3)
    χ = B[x, :]

    figure("quadrature points", figsize=(7,9))
    clf()
    subplot(411)
    plot(x, χ, length(x) < 30 ? ".-" : "-")
    margins(0.1, 0.1)
    no_tick_labels()
    subplot(412)
    rplot(t)
    margins(0.1, 0.1)
    xl = xlim()
    no_tick_labels()
    ylabel("Multiplicity")
    subplot(413)
    plot(B.x, 1:length(B.x), ".-")
    margins(0.1,0.1)
    xlim(xl)
    no_tick_labels()
    ylabel(L"i")
    subplot(414)
    plot(B.x, B.w, ".-")
    margins(0.1,0.1)
    xlim(xl)
    xlabel(L"x")
    ylabel(L"w_i")
    tight_layout()
    savefig("docs/src/figures/quadrature-points.svg")
end

function function_interpolation()
    t = LinearKnotSet(7, 0, 7, 10)
    B = BSpline(t)
    x = range(0, stop=7, length=301)
    χ = B[x, :]
    c = B \ sin
    i = 1:size(B,2)

    figure("function interpolation",figsize=(7,9))
    clf()
    subplot(311)
    plot(x, χ*c)
    plot([mean_position(x, view(χ, :, j)) for j in i], c, "s-")
    no_tick_labels()
    subplot(312)
    plot(x, χ*c-sin.(x))
    ylabel("Error")
    no_tick_labels()
    subplot(313)
    plot(x, χ)
    xlabel(L"x")
    tight_layout()
    savefig("docs/src/figures/function-interpolation.svg")
end

function restricted_basis_interpolation()
    t = LinearKnotSet(7, 0.0, 1.0, 6)
    x = range(first(t), stop=last(t), length=300)[2:end-1]

    B = BSpline(t,3)
    B̃ = B[:,2:end-1]

    f1 = x -> sin(2π*x)
    f2 = x -> cos(2π*x)

    c1 = B \ f1
    c2 = B \ f2

    c̃1 = B̃ \ f1
    c̃2 = B̃ \ f2

    χ = B[x, :]
    χ̃ = B̃[x, :]
    x_avg = [mean_position(x, view(χ, :, j)) for j in 1:size(B,2)]
    x̃_avg = [mean_position(x, view(χ̃, :, j)) for j in 1:size(B̃,2)]

    figure("restricted basis interpolation", figsize=(7,9))
    clf()
    subplot(321)
    l=plot(x, χ*c1)[1]
    plot(x_avg, c1, "s:", color=l.get_color())
    l=plot(x, χ*c2)[1]
    plot(x_avg, c2, "s:", color=l.get_color())
    no_tick_labels()
    subplot(322)
    l=plot(x, χ̃*c̃1)[1]
    plot(x̃_avg, c̃1, "s:", color=l.get_color())
    l=plot(x, χ̃*c̃2)[1]
    plot(x̃_avg, c̃2, "s:", color=l.get_color())
    axes_labels_opposite(:y)
    no_tick_labels()
    subplot(324)
    semilogy(x, abs.(χ̃*c̃1 - f1.(x)))
    semilogy(x, abs.(χ̃*c̃2 - f2.(x)))
    axes_labels_opposite(:y)
    yl = ylim()
    no_tick_labels()
    ylabel("Error")
    subplot(323)
    semilogy(x, abs.(χ*c1 - f1.(x)))
    semilogy(x, abs.(χ*c2 - f2.(x)))
    no_tick_labels()
    ylabel("Error")
    ylim(yl)
    subplot(325)
    plot(x, χ)
    xlabel(L"x")
    yl = ylim()
    subplot(326)
    plot(x, χ̃)
    axes_labels_opposite(:y)
    ylim(yl)
    xlabel(L"x")
    tight_layout()
    savefig("docs/src/figures/restricted-basis-interpolation.svg")
end

function smooth_interpolation()
    f = x -> sin(2π*x)

    rng = MersenneTwister(123);

    N = 10
    x = clamp.(sort(range(0, stop=1, length=N) + 0.1(2rand(rng,N) .- 1)), 0, 1);
    y = f.(x) + 0.1(2rand(rng,N) .- 1);

    t3 = LinearKnotSet(3, 0.0, 1.0, 6);
    t4 = LinearKnotSet(4, 0.0, 1.0, 6);
    B3 = BSpline(t3,0)
    B4 = BSpline(t4,0)

    c3 = B3[x,:] \ y
    c4 = B4[x,:] \ y

    r = range(first(t3), stop=last(t3), length=300)
    χ3 = B3[r,:]
    χ4 = B4[r,:]

    figure("smooth interpolation")
    clf()
    plot(x, y, "s", label="Samples")
    plot(r, χ3*c3, label="3rd order spline")
    plot(r, χ4*c4, label="4th order spline")
    plot(r, f.(r), "--", label=L"\sin(2\pi x)")
    legend()
    tight_layout()
    savefig("docs/src/figures/smooth-interpolation.svg")
end

function diagonal_operators()
    k = 7
    N = 31

    a,b = 0,70

    coulomb(r) = -1/r

    cfigure("V(x)", figsize=(7,9)) do
        for (j,(t,x)) in enumerate([(LinearKnotSet(k, a, b, N),
                                     range(a, stop=b, length=500)[2:end-1]),
                                    (ExpKnotSet(k, -1.0, log10(b), N),
                                     10 .^ range(-1.0, stop=log10(b), length=500)[2:end-1])])
            B = BSpline(t,3)[:,2:end-1]
            S = B'B

            χ = B[x,:]

            f = B \ x -> x^2*exp(-x)
            g = B \ x -> -x*exp(-x)

            V = Matrix(coulomb, B)
            g̃ = S \ V*f

            csubplot(3,2,(j-1)+1, nox=true) do
                plot(x, χ*f, label=L"f(x)")
                plot(x, χ*g̃, label=L"\tilde{g}(x)")
                plot(x, χ*g, "--", label=L"g(x)")
                yl=ylim()
                plot(x, coulomb.(x), label=L"V(x)")
                ylim(yl)
                legend(framealpha=0.75)
                xscale("log")
                iseven(j) && axes_labels_opposite(:y)
            end
            csubplot(3,2,(j-1)+3, nox=true) do
                plot(x, χ*(g-g̃), label=L"g(x)-\tilde{g}(x)")
                legend(framealpha=0.75)
                xscale("log")
                ylabel("Error")
                iseven(j) && axes_labels_opposite(:y)
            end
            csubplot(3,2,(j-1)+5) do
                plot(x, χ)
                xscale("log")
                xlabel(L"x")
                iseven(j) && axes_labels_opposite(:y)
            end
        end
    end
    savefig("docs/src/figures/diagonal-operators.svg")
end

function find_second_derivative(B, f::Function)
    S = B'B
    D = Derivative(axes(B,1))
    ∇² = B'D'D*B

    # Project function onto B-spline basis
    cf = B \ f
    # Find derivative
    cg = S \ ∇²*cf

    cf,cg
end

function sine_derivative()
    t = LinearKnotSet(10, 0, 10, 30);
    B = BSpline(t)[:,2:end-1]

    f = x -> sin(2π*x)
    g = x -> -4π^2*sin(2π*x)

    cf,cg = find_second_derivative(B, f)

    x = range(first(t), stop=last(t), length=1001)
    χ = B[x, :]

    cfigure("derivatives", figsize=(7,9)) do
        csubplot(411,nox=true) do
            l=plot(x, χ*cf, label=L"\tilde{f}(x)")[1]
            plot(x, f.(x), ":", label=L"f(x)")
            legend(framealpha=0.75)
        end
        csubplot(412,nox=true) do
            l=plot(x, (χ*cg), label=L"\tilde{g}(x)")[1]
            plot(x, g.(x), ":", label=L"g(x)")
            legend(framealpha=0.75)
        end
        csubplot(413,nox=true) do
            semilogy(x, abs.(χ*cg-g.(x)), label=L"|\tilde{g}(x)-g(x)|")
            legend(framealpha=0.75)
            ylim(1e-5,1)
            ylabel("Error")
        end
        csubplot(414) do
            plot(x, χ)
            xlabel(L"x")
        end
    end
    savefig("docs/src/figures/sine-derivative.svg")
end

function ode_hookes_law(xₘₐₓ, kspring, k, N)
    t = LinearKnotSet(k, 0, xₘₐₓ, N)
    # By omitting the first basis function, we enforce V(0) = 0
    B = BSpline(t,0)[:,2:end]
    S = B'B

    D = Derivative(axes(B, 1))
    ∇ = B'D*B

    # Hooke's law
    F = x -> -kspring*x
    # Exact potential
    V = x -> kspring*x^2/2

    # Expand Hooke's law on B-splines
    cF = B \ F
    # Solve for expansion coefficients of potential
    cV = -∇ \ S*cF

    x = range(first(t), stop=last(t), length=500)
    χ = B[x,:]
    x_avg = [mean_position(x, view(χ, :, j)) for j in 1:size(B,2)]

    cfigure("Hooke's law",figsize=(7,9)) do
        csubplot(411,nox=true) do
            l=plot(x, χ*cF, label=L"\tilde{F}(x)")[1]
            plot(-x, -χ*cF, "--", color=l.get_color())
            plot(x_avg, cF, ".:", color=l.get_color(), label=L"c_F")
            plot(x, F.(x), ":", label=L"F(x)")
            legend(framealpha=0.75)
        end
        csubplot(412,nox=true) do
            l=plot(x, χ*cV, label=L"\tilde{V}(x)")[1]
            plot(-x, χ*cV, "--", color=l.get_color())
            plot(x_avg, cV, ".:", color=l.get_color(), label=L"c_V")
            plot(x, V.(x), ":", label=L"V(x)")
            legend(framealpha=0.75)
        end
        csubplot(413,nox=true) do
            plot(x, χ*cV - V.(x))
            xl = xlim()
            xlim(-xl[2],xl[2])
            ylabel("Error")
        end
        csubplot(414) do
            plot(x, χ)
            xl = xlim()
            xlim(-xl[2],xl[2])
            xlabel(L"x")
        end
    end
    savefig("docs/src/figures/hookes-law-$(k)-$(N).svg")
end

struct ShiftAndInvert{TA,TB,TT}
    A⁻¹::TA
    B::TB
    temp::TT
end

Base.size(S::ShiftAndInvert, args...) = size(S.A⁻¹, args...)
Base.eltype(S::ShiftAndInvert) = eltype(S.A⁻¹)

function LinearAlgebra.mul!(y,M::ShiftAndInvert,x)
    mul!(M.temp, M.B, x)
    ldiv!(y, M.A⁻¹, M.temp)
end

construct_linear_map(A,B,σ=0) =
    ShiftAndInvert(factorize(A-σ*B),B,Vector{eltype(A)}(undef, size(A,1)))

function hydrogen_eigenstates()

    k = 7
    N = 31

    a,b = 0,70

    coulomb(r) = -1/r

    nev = 5
    σ = -0.5

    n = 1:nev

    cfigure("Hydrogen", figsize=(7,9)) do
        for (j,(t,x,tol)) in enumerate([(LinearKnotSet(k, a, b, N),
                                         range(a, stop=b, length=500)[2:end-1],
                                         9e-3),
                                        (ExpKnotSet(k, -1.0, log10(b), N),
                                         10 .^ range(-1.0, stop=log10(b), length=500)[2:end-1],
                                         2e-7)])
            B = BSpline(t,3)[:,2:end-1]
            S = B'B

            χ = B[x,:]

            D = Derivative(axes(B, 1))

            ∇² = B'D'D*B

            T = -∇²/2

            V = Matrix(coulomb, B)

            H = T + V

            schurQR,history = partialschur(construct_linear_map(H, S, σ), nev=nev)
            println(history)

            θ = schurQR.eigenvalues
            E = real(σ .+ inv.(θ))

            csubplot(3,2,(j-1)+1, nox=true) do
                for i = 1:nev
                    ϕ = schurQR.Q[:,i]
                    plot(x, real(E[i]) .+ abs2.(χ*ϕ)/(ϕ'S*ϕ)[1])
                end
                # legend(framealpha=0.75)
                xscale("log")
                iseven(j) && axes_labels_opposite(:y)
            end
            csubplot(3,2,(j-1)+3, nox=true) do
                plot(x, coulomb.(x))
                xscale("log")
                iseven(j) && axes_labels_opposite(:y)
            end
            csubplot(3,2,(j-1)+5) do
                plot(x, χ)
                plot(x, χ)
                xscale("log")
                xlabel(L"x")
                iseven(j) && axes_labels_opposite(:y)
            end
        end
    end
    savefig("docs/src/figures/hydrogen-eigenstates.svg")
end

mkpath("docs/src/figures")
cardinal_splines()
discontinuous_splines()
full_multiplicity_splines()
spline1d()
spline2d()
logo()
quadrature_points()
function_interpolation()
restricted_basis_interpolation()
smooth_interpolation()
diagonal_operators()
sine_derivative()
ode_hookes_law(3, 0.1, 7, 30)
ode_hookes_law(3, 0.1, 3, 1)
hydrogen_eigenstates()
