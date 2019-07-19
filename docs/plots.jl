using PyPlot
using Jagot.plotting
plot_style("ggplot")
using PyPlotRecipes
using PyCall
PathEffects = pyimport("matplotlib.patheffects")
using Statistics
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
    figure("splines",figsize=(7,9))
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

mkpath("docs/src/figures")
cardinal_splines()
discontinuous_splines()
full_multiplicity_splines()
spline1d()
spline2d()
logo()
