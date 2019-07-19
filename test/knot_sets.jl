using IntervalSets
import BSplinesQuasi: RightContinuous, within_interval, within_support

"""
    within_interval_linear(x, interval)

Return the range of indices of elements of the sorted vector/range `x`
that fall within `interval`. This loops through all elements of `x`
(linear complexity), so it should not be used for large number of
elements.
"""
function within_interval_linear(x, interval)
    @assert issorted(x)
    N = length(x)+1
    i = N
    j = 0
    for (k,e) in enumerate(x)
        if e ∈ interval
            if k < i
                i = k
            else
                j = k
            end
        end
    end

    i:j
end

function test_within_interval(x, interval, expected=nothing)
    linear_expected = within_interval_linear(x, interval)
    if isnothing(expected)
        expected = linear_expected
    else
        # This is mostly a sanity check to test if the expectations
        # are actually correct.
        expected ≠ linear_expected &&
            throw(ArgumentError("Invalid expectations: $(expected), actual: $(linear_expected)"))
    end

    result = :(within_interval($x, $interval))
    expr = :($result == $expected)
    if !(@eval $expr)
        println("Looking for elements of $x ∈ $interval, got $(@eval $result), expected $expected")
        length(x) < 30 && println("    x = ", collect(enumerate(x)), "\n")
    end
    @eval @test $expr
end

@testset "Knot sets" begin
    @testset "k = $k" for k ∈ 1:6
        @testset "Full multiplicity" begin
            t = LinearKnotSet(k, 0, 1, 1)
            @test length(t) == 2k
            @test collect(t) == vcat(fill(0, k), fill(1, k))
        end
        @testset "Simple endpoints" begin
            t = LinearKnotSet(k, 0, 1, k, 1, 1)
            @test length(t) == k+1
            @test collect(t) == range(0, stop=1, length=k+1)
        end
    end

    @testset "Invalid order/multiplicities" begin
        @testset "Invalid order" begin
            @test_throws ArgumentError LinearKnotSet(0, 0, 1, 4)
            @test_throws ArgumentError LinearKnotSet(-1, 0, 1, 4)
        end
        @testset "Invalid multipicities" begin
            @test_throws ArgumentError LinearKnotSet(1, 0, 1, 4, 0, 1)
            @test_throws ArgumentError LinearKnotSet(1, 0, 1, 4, 1, 0)
            @test_throws ArgumentError LinearKnotSet(1, 0, 1, 4, 1, 2)
        end
        @testset "Not enough interior points" begin
            @test_throws ArgumentError LinearKnotSet(2, 0, 1, 1, 1, 1)
        end
    end

    @testset "Properties" begin
        @testset "#intervals = $N" for N = 1:4
            @testset "k = $k" for k = 1:N
                @testset "Simple multiplicities" begin
                    t = LinearKnotSet(k, 0.0, 1.0, N, 1, 1)
                    @test numintervals(t) == N
                    @test numfunctions(t) == N - k + 1
                end
            end
        end
    end

    @testset "Compact support" begin
        a,b = 0,1
        x = range(a, stop=b, length=21)
        @testset "Interval coverage" begin
            @testset "Two intervals" begin
                test_within_interval(x, 0..0.5, 1:11)
                test_within_interval(x, RightContinuous(0,0.5))
                test_within_interval(x, RightContinuous(0,0.5))
                test_within_interval(x, RightContinuous(0.25,0.5))
                test_within_interval(x, RightContinuous(0.25,0.5))
            end
            @testset "Three intervals" begin
                test_within_interval(x, RightContinuous(0,1/3), 1:7)
                test_within_interval(x, RightContinuous(1/3,2/3), 8:14)
                test_within_interval(x, 2/3..1, 15:21)
            end
            @testset "Open interval" begin
                test_within_interval(x, OpenInterval(0.2,0.4), 6:8)
            end
            @testset "Random intervals" begin
                @testset "L=$L" for L=[:closed,:open]
                    @testset "R=$R" for R=[:closed,:open]
                        for i = 1:1 # 20
                            interval = Interval{L,R}(minmax(rand(),rand())...)
                            test_within_interval(x, interval)
                        end
                    end
                end
            end
        end

        @testset "Support of Heavyside splines" begin
            t = LinearKnotSet(1, a, b, 4)
            supports = [within_support(x, t, j)
                        for j = 1:numfunctions(t)]
            @test length(supports) == 4
            # Each basis function should cover one interval only (since order = 1).
            @test all(length.(supports) .== 1)
            # Test that all elements of x are covered by the basis
            # functions, and that none of the basis functions overlap.
            @test first.(first.(supports)) == [1:5, 6:10, 11:15, 16:21]
        end

        @testset "Support of linear splines" begin
            @testset "Simple multiplicity" begin
                t = LinearKnotSet(2, a, b, 2, 1, 1)
                supports = [within_support(x, t, j)
                            for j = 1:numfunctions(t)]
                @test length(supports) == 1
                @test supports[1] == [(1:10,1), (11:21,2)]
            end
            @testset "Full multiplicity" begin
                @testset "#intervals = 2" begin
                    t = LinearKnotSet(2, a, b, 2)
                    supports = [within_support(x, t, j)
                                for j = 1:numfunctions(t)]
                    @test length(supports) == 3
                    @test supports[1] == [(1:10,2)]
                    @test supports[2] == [(1:10,2), (11:21,3)]
                    @test supports[3] == [(11:21,3)]
                end
                @testset "#intervals = 3" begin
                    t = LinearKnotSet(2, a, b, 3)
                    supports = [within_support(x, t, j)
                                for j = 1:numfunctions(t)]
                    @test length(supports) == 4
                    @test supports[1] == [(1:7,2)]
                    @test supports[2] == [(1:7,2), (8:14,3)]
                    @test supports[3] == [(8:14,3), (15:21,4)]
                    @test supports[4] == [(15:21,4)]
                end
            end
        end
    end
end
