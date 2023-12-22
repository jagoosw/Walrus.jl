@testset "Wind stress" begin
    grid = RectilinearGrid(size = (2, 2, 2), extent = (2, 2, 2))

    wind_stress_boundary_conditions = WindStressBoundaryConditions(; reference_wind_speed = 1., reference_wind_direction = 90.)

    model = NonhydrostaticModel(; grid, boundary_conditions = (u = FieldBoundaryConditions(top = wind_stress_boundary_conditions.u),
                                                               v = FieldBoundaryConditions(top = wind_stress_boundary_conditions.v)))

    for n=1:100000
        time_step!(model, 1)
    end

    @test -1 <= model.velocities.u[1, 1, 2] < 0.1
    @test model.velocities.u[1, 1, 1] ≈ 0
    @test model.velocities.v[1, 1, 2] ≈ 0
end