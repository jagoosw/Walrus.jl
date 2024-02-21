@testset "Wind stress" begin
    grid = RectilinearGrid(arch; size = (2, 2, 2), extent = (2, 2, 2))

    wind_stress_boundary_conditions = WindStressBoundaryConditions(; reference_wind_speed = 0., reference_wind_direction = 90.)

    model = NonhydrostaticModel(; grid, boundary_conditions = (u = FieldBoundaryConditions(top = wind_stress_boundary_conditions.u),
                                                               v = FieldBoundaryConditions(top = wind_stress_boundary_conditions.v)))

    for n=1:1000
        time_step!(model, 1)
    end

    @test all(Array(model.velocities.u) .≈ 0) # no wind no stress

    wind_stress_boundary_conditions = WindStressBoundaryConditions(; reference_wind_speed = 1., reference_wind_direction = 90.)

    model = NonhydrostaticModel(; grid, boundary_conditions = (u = FieldBoundaryConditions(top = wind_stress_boundary_conditions.u),
                                                               v = FieldBoundaryConditions(top = wind_stress_boundary_conditions.v)))

    time_step!(model, 1)
    
    for n=1:10000
        time_step!(model, 1)
    end

    @test -1 <= Array(model.velocities.u, 1, 1, 2) < 0.1
    @test Array(interior(model.velocities.u, 1, 1, 1))[1] ≈ 0
    @test all(Array(model.velocities.v) .≈ 0)

    # precomputed roughness lengths

    wind_stress_boundary_conditions = WindStressBoundaryConditions(; reference_wind_speed = 0., 
                                                                     reference_wind_direction = 90., 
                                                                     drag_coefficient = 
                                                                        LogarithmicNeutralWind(; precomputed_roughness_length = true, arch))

    model2 = NonhydrostaticModel(; grid, boundary_conditions = (u = FieldBoundaryConditions(top = wind_stress_boundary_conditions.u),
                                                                v = FieldBoundaryConditions(top = wind_stress_boundary_conditions.v)))

    for n=1:1000
        time_step!(model2, 1)
    end

    @test all(Array(model2.velocities.u) .≈ 0) # no wind no stress

    wind_stress_boundary_conditions = WindStressBoundaryConditions(; reference_wind_speed = 1., 
                                                                     reference_wind_direction = 90., 
                                                                     drag_coefficient = 
                                                                        LogarithmicNeutralWind(; precomputed_roughness_length = true, arch))

    model2 = NonhydrostaticModel(; grid, boundary_conditions = (u = FieldBoundaryConditions(top = wind_stress_boundary_conditions.u),
                                                                v = FieldBoundaryConditions(top = wind_stress_boundary_conditions.v)))

    time_step!(model2, 1)

    for n=1:10000
        time_step!(model2, 1)
    end

    @test -1 <= Array(interior(model2.velocities.u, 1, 1, 2))[1] < 0.1
    @test Array(interior(model2.velocities.u, 1, 1, 1))[1] ≈ 0
    @test all(Array(model2.velocities.v) .≈ 0)

    @test all(Array(model.velocities.u) .≈ Array(model2.velocities.u))
end