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

function Bspline_text(x, ϕ, j, k, color)
    xⱼ = ϕ'*Diagonal(x)*ϕ/(ϕ'ϕ)
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

mkpath("docs/src/figures")
cardinal_splines()
discontinuous_splines()
full_multiplicity_splines()
