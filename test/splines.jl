@testset "Splines" begin
    @testset "Discontinuous knot set" begin
        t = ArbitraryKnotSet(3, [0.0, 1, 1, 3, 4, 6], 1, 3)
        B = BSpline(t)

        @test all(B[0.0, :] .== 0)
        @test all(B[0.0:0.5:1.0, :][1, :] .== 0)
    end

    @testset "Evaluate B-splines" begin
        k = 3
        t = LinearKnotSet(k, 0, 1, 2)
        B = BSpline(t)
        @testset "Eval on subintervals" begin
            x₁ = range(-1,stop=-0.5,length=10)
            χ₁ = B[x₁, :]
            @test norm(χ₁) == 0

            function testbasis(x)
                χ = B[x, :]
                B̃ = spzeros(Float64, length(x), 4)
                B̃[:,1] = (x .>= 0) .* (x .< 0.5) .* ((2x .- 1).^2)
                B̃[:,2] = (x .>= 0) .* (x .< 0.5) .* (2/3*(1 .- (3x .- 1).^2)) +
                    (x .>= 0.5) .* (x .< 1)  .* (2*(x .- 1).^2)
                B̃[:,3] = (x .>= 0) .* (x .< 0.5) .* (2*x.^2) +
                    (x .>= 0.5) .* (x .< 1)  .* (2/3*(1 .- (3x .- 2).^2))
                B̃[:,4] = (x .>= 0.5) .* (x .<= 1) .* ((2x .- 1).^2)
                for j = 1:4
                    δ,δr = vecdist(χ[:,j],B̃[:,j])
                    @test δ < 1e-15
                    @test δr < 1e-15
                end
            end
            testbasis(range(0,stop=1,length=50))
            testbasis(range(-0.5,stop=0.6,length=40))
            testbasis(range(0.5,stop=1.6,length=40))
        end
    end

    @testset "Evaluate spline" begin
        @testset "Linear knot set" begin
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
        @testset "Discontinuous knot set" begin
            t = ArbitraryKnotSet(3, [0.0, 1, 1, 3, 4, 6], 1, 3)
            B = BSpline(t,3)

            χ = B[B.x, :]
            @test all(isfinite, χ)
            @test all(!isnan, χ)
        end
    end
end
