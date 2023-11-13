using Resparse
using Test
using  ForwardDiff,GenericLinearAlgebra,LinearAlgebra

@testset "Resparse.jl" begin
    # Write your tests here.

    A = randn(2, 2)
    f(x) = svd(x).S |> sum
    g(x) = Resparse.svd_lowrank(x, 2, 2).S |> sum
    @time ForwardDiff.gradient(g, A)

end
