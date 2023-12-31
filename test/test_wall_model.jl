@testset "Wall model" begin
    grid = RectilinearGrid(size = (2, 2, 2), extent = (2, 2, 2))

    stress_boundary_conditions = WallStressBoundaryConditions()

    model = NonhydrostaticModel(; grid, boundary_conditions = (u = FieldBoundaryConditions(bottom = stress_boundary_conditions.u),
                                                               v = FieldBoundaryConditions(bottom = stress_boundary_conditions.v)))

    for n=1:100
        time_step!(model, 1)
    end

    @test all(model.velocities.u .≈ 0) && all(model.velocities.v .≈ 0) && all(model.velocities.w .≈ 0) # no velocity change when no velocity

    set!(model, u = 1)

    for n=1:1000
        time_step!(model, 1)
    end

    @test all(interior(model.velocities.u, :, :, 1) .< 1) # when moving, it is slowed
end