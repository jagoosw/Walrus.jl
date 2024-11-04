using Walrus.Interpolations: Limited

@testset "Interpolaiton" begin
    x = [0:10^-6:2π;]
    y = cos.(x)

    itp = SimpleInterpolation(x, y; boundary_condition = Limited())

    test_x = rand(0:10^-6:2π, 1000)

    # atol is the Rolle theory maxiumum error for Limited interpolation + a little
    @test all([isapprox(itp(x), cos(x), atol = 2 * abs(0.01^2 / 8 * cos(floor(x; digits = 2)))) for x in test_x])
    
    test_x = rand(0:10^-6:4π, 1000)

    itp = SimpleInterpolation(x, y)

    @test all([isapprox(itp(x), cos(x), atol = 15 * abs(0.01^2 / 8 * cos(floor(x; digits = 2)))) for x in test_x]) 
end