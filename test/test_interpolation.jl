@testset "Interpolaiton" begin
    Δx = 1e-4
    x = 0:Δx:1 - Δx
    y = @. sin(2π * x)

    itp = SimpleInterpolation(Array(x), y)

    test_x = rand(1000) .* 2

    # atol is the Rolle theory maxiumum error for linear interpolation
    # https://en.wikipedia.org/wiki/Linear_interpolation#Linear_interpolation_as_an_approximation
    abs_f′′_x₀_x₁ = @. (2π)^2 * max(abs(sin(2π * ceil(test_x; digits = ceil(Int, -log10(Δx))))), abs(sin(2π * floor(test_x; digits = ceil(Int, -log10(Δx))))))

    @test all([isapprox(itp(x), sin(2π * x), atol = Δx^2 / 8 * abs_f′′_x₀_x₁[n]) for (n, x) in enumerate(test_x)]) 
end