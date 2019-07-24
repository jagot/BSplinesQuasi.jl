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

function mean_color(color::String)
    c = parse(Colorant, color)
    mean([c.r,c.g,c.b])
end

lerp(a,b,t) = (1-t)*a + t*b

mean_position(x, ϕ) = ϕ'*Diagonal(x)*ϕ/(ϕ'ϕ)

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
