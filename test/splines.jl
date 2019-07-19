@testset "Splines" begin
    @testset "Discontinuous knot set" begin
        t = ArbitraryKnotSet(3, [0.0, 1, 1, 3, 4, 6], 1, 3)
        B = BSpline(t)

        @test all(B[0.0, :] .== 0)
        @test all(B[0.0:0.5:1.0, :][1, :] .== 0)
    end

    @testset "Evaluate spline" begin
        k = 4
        t = LinearKnotSet(k, 0, 1, 5)
        B = BSpline(t)
        x = range(0, stop=1, length=301)
        χ = B[x, :]

        i = 1:size(B,2)
        c = sin.(i)

        # Whether we evaluate the spline or the basis functions should
        # not make a difference.
        s = (B*c)[x]
        s′ = χ*c
        @test all(isfinite, s)
        @test all(isfinite, s′)
        @test all(!isnan, s)
        @test all(!isnan, s′)
        @test s == s′
    end
end
