@testset "Operators" begin
    @testset "Coulomb operator" begin
        k = 7
        N = 31

        a,b = 0,70

        coulomb(r) = -1/r
        @testset "$name knot set" for (name, t, tol) in [("Linear", LinearKnotSet(k, a, b, N), 9e-3),
                                                         ("Exponential", ExpKnotSet(k, -1.0, log10(b), N), 2e-7)]
            B = BSpline(t,3)[:,2:end-1]
            S = B'B

            f = B \ x -> x^2*exp(-x)
            g = B \ x -> -x*exp(-x)

            V = Matrix(coulomb, B)
            g̃ = S \ V*f

            @test g ≈ g̃ atol=tol
        end
    end
end
