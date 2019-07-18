import Base: first, last, length,
    getindex, lastindex, eachindex, iterate,
    eltype, similar,
    show

# * AbstractKnotSet
"""
    AbstractKnotSet{T,k,ml,mr}

Abstract base for B-spline knot sets. `T` is the `eltype` of the knot
set, `k` is the order of the piecewise polynomials (order = degree +
1) and `ml` and `mr` are the knot multiplicities of the left and right
endpoints, respectively.
"""
abstract type AbstractKnotSet{k,ml,mr,T} <: AbstractVector{T} end

function assert_multiplicities(k,ml,mr,t)
    k > 0 || throw(ArgumentError("Polynomial order must be a positive integer, got k = $k"))
    ml ∈ 1..k || throw(ArgumentError("Left multiplicity has to be in range 1..$(k), got $ml"))
    mr ∈ 1..k || throw(ArgumentError("Right multiplicity has to be in range 1..$(k), got $mr"))
    n = length(t)
    nreq = max(2, k + 3 - ml - mr)
    n ≥ nreq || throw(ArgumentError("Polynomial order k = $(k) and left,right multiplicities $(ml),$(mr) require a knot sequence of at least $(nreq) points"))
end

order(t::AbstractKnotSet{k}) where k = k
numintervals(t::AbstractKnotSet) = length(t.t)-1
numfunctions(t::AbstractKnotSet) = length(t) - order(t)

leftmultiplicity(::AbstractKnotSet{k,ml,mr}) where {k,ml,mr} = ml
rightmultiplicity(::AbstractKnotSet{k,ml,mr}) where {k,ml,mr} = mr

first(t::AbstractKnotSet) = first(t.t)
last(t::AbstractKnotSet) = last(t.t)
length(t::AbstractKnotSet) = length(t.t) + leftmultiplicity(t) + rightmultiplicity(t) - 2
size(t::AbstractKnotSet) = (length(t),)

function getindex(t::AbstractKnotSet{k,ml,mr}, i::Integer) where {k,ml,mr}
    i < 1 && throw(BoundsError("Trying to access knot set of length $(length(t)) at index $i"))

    ni = numintervals(t)
    if i < ml
        first(t)
    elseif i < ml + ni
        t.t[i-ml+1]
    elseif i < ml + ni + mr
        last(t)
    else
        throw(BoundsError("Trying to access knot set of length $(length(t)) at index $i"))
    end
end
lastindex(t::AbstractKnotSet) = length(t)
eachindex(t::AbstractKnotSet) = 1:length(t)
# function iterate(t::AbstractKnotSet, (element,i)=(t[1],0))
#     i ≥ length(t) && return nothing
#     element, (t[i+1], i+1)
# end

eltype(t::AbstractKnotSet{k,ml,mr,T}) where {k,ml,mr,T} = T
similar(t::AbstractKnotSet{k,ml,mr,T}) where {k,ml,mr,T} = Vector{T}(undef, length(t))

function show_order_multiplicities(io::IO, t::AbstractKnotSet{k,ml,mr}) where {k,ml,mr}
    write(io, "order k = $(k)")
    k < 7 && write(io, " (",
                  ["constant", "linear", "parabolic", "cubic", "quartic", "quintic"][k],
                  ")")
    ml == k && mr == k && return
    write(io, " (left, right multiplicities: $(ml), $(mr))")
end

function show(io::IO, t::AbstractKnotSet)
    write(io, "$(typeof(t).name)($(eltype(t))) of ")
    show_order_multiplicities(io, t)
    write(io, " on $(first(t)..last(t)) ($(numintervals(t)) intervals)")
end

function find_interval(t::AbstractKnotSet{T,k,ml,mr}, x, i=ml) where {T,k,ml,mr}
    (x < first(t) || x > last(t)) && return nothing
    x == last(t) && return ml + length(t.t) - 1
    for r ∈ i:length(t)
        t[r] > x && return r-1
    end
    @assert false
end

const RightContinuous{T} = Interval{:closed,:open,T}

function lfloor(::Type{T}, r) where T
    i = floor(T, r)
    i == r ? i - 1 : i
end

"""
    within_interval(x, interval)

Return the indices of the elements of `x` that lie within the given
closed `interval`.
"""
function within_interval(x::AbstractRange, interval::Interval{L,R}) where {L,R}
    N = length(x)
    δx = step(x)
    l = leftendpoint(interval)
    r = rightendpoint(interval)
    a = max(1, ceil(Int, (l-x[1])/δx) + 1)
    b = min(N, ceil(Int, (r-x[1])/δx))
    a + (x[a] == l && L == :open):b + (b < N && x[b+1] == r && R == :closed)
end

"""
    within_support(x, t, j)

Return the indices of the elements of `x` that lie withing the compact
support of the `j`th basis function (enumerated `1..n`), given the
knot set `t`. For each index of `x` that is covered, the index `k` of
the interval within which `x[i]` falls is also returned.
"""
function within_support(x::AbstractRange, t::AbstractKnotSet, j::Integer)
    isempty(x) && return 1:0
    k = order(t)
    # The last interval includes the right endpoint as well, whereas
    # the other intervals are only right-continuous.
    IntervalKind(i) = i == numintervals(t) ? ClosedInterval : RightContinuous
    ml = leftmultiplicity(t)
    # @show j, j:j+k-1
    supports = [(within_interval(x, IntervalKind(i-ml+1)(t[i], t[i+1])), i)
                for i = j:j+k-1]
    filter(s -> !isempty(s[1]), supports)
end

# * Specific knot sets
# ** Arbitrary knot set
struct ArbitraryKnotSet{k,ml,mr,T,V<:AbstractVector{T}} <: AbstractKnotSet{k,ml,mr,T}
    t::V
    ArbitraryKnotSet{k,ml,mr}(t::V) where {k,ml,mr,T,V<:AbstractVector{T}} =
        new{k,ml,mr,T,V}(sort(t))
end
function ArbitraryKnotSet(k::Integer, t::AbstractVector,
                          ml::Integer=k, mr::Integer=k)
    assert_multiplicities(k, ml, mr, t)
    ArbitraryKnotSet{k,ml,mr}(t)
end

# ** Linear knot set
struct LinearKnotSet{k,ml,mr,T,R<:AbstractRange{T}} <: AbstractKnotSet{k,ml,mr,T}
    t::R
    LinearKnotSet{k,ml,mr}(t::R) where {k,ml,mr,T,R<:AbstractRange{T}} =
        new{k,ml,mr,T,R}(t)
end
function LinearKnotSet(k::Integer, a, b, N::Integer,
                       ml::Integer=k, mr::Integer=k)
    t = range(a, stop=b, length=N+1)
    assert_multiplicities(k, ml, mr, t)
    LinearKnotSet{k,ml,mr}(t)
end

# ** Exponential knot set
struct ExpKnotSet{k,ml,mr,T,R<:AbstractRange{T},TV<:AbstractVector{T}} <: AbstractKnotSet{k,ml,mr,T}
    exponents::R
    base::T
    t::TV
    include0::Bool
end
function ExpKnotSet(k::Integer, a::T, b::T, N::Integer,
                    ml::Integer=k, mr::Integer=k;
                    base::T=T(10), include0::Bool=true) where T
    exponents = range(a, stop=b, length=include0 ? N : N+1)
    t = base .^ exponents
    assert_multiplicities(k, ml, mr, t)
    ExpKnotSet{k,ml,mr}(exponents, eltype(t)(base), include0 ? vcat(0,t) : t, include0)
end

function show(io::IO, t::ExpKnotSet)
    write(io, "$(typeof(t).name)($(eltype(t))) of  on ")
    show_order_multiplicities(io, t)
    write(io, " on ")
    t.include0 && write(io, "0,")
    write(io, "$(t.base^first(t.exponents))..$(t.base^last(t.exponents)) ($(numintervals(t)) intervals)")
end

# function arcsin_knot_set(k::Integer, a::Integer, b::Integer, N::Integer)
#     N2 = N/2
#     arcsin.(range(-N2, stop=N2, lengthN+1)/N2)*(b-a)/π+(a+b)/2
# end

# arcsin_half_knot_set(a, b, N) = arcsin(linspace(0,N,N+1)/N)*(b-a)*2.0/π+a

# function exp_linear_knot_set(a, b, N)
#     ap = a != 0 ? a : 1e-1
#     @assert t[2][1]+t[2][2] == N
#     [logspace(log10(ap),log10(t[1]), t[2][1]); linspace(t[1],b,t[2][2]+1)[2:end]]
# end

@recipe function plot(t::AbstractKnotSet)
    markershape --> :circle
    y = similar(t)
    y′ = 1
    t′ = -Inf
    for (i,e) in enumerate(t)
        y′ = e == t′ ? (y′ + 1) : 1
        t′ = e
        y[i] = y′
    end
    collect(t),y
end

export ArbitraryKnotSet, LinearKnotSet, ExpKnotSet, order, numintervals, numfunctions
