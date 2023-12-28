@testset "Interpolaiton" begin
    x = [0:0.01:1;]
    y = cos.(x)

    itp = SimpleInterpolation(x, y)

    test_x = rand(1000)

    # atol is the Rolle theory maxiumum error for linear interpolation
    @test all([isapprox(itp(x), cos(x), atol = 0.01^2 / 8 * cos(floor(x; digits = 2))) for x in test_x]) 
end