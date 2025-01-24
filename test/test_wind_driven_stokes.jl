using Oceananigans.Units

@testset "Wind stress" begin
    grid = RectilinearGrid(arch; size = (32, 32, 16), extent = (6.55*3, 6.55*3, 10))

    wind_stress_boundary_conditions = WindStressBoundaryConditions(; reference_wind_speed = 1., reference_wind_direction = 270., drag_coefficient = LogarithmicNeutralWind(; precompute_drag_coefficients = true))
    stokes_drift = WindDrivenStokesDriftSetup(; wind = wind_stress_boundary_conditions.u.condition.func,
                                                direction = 270,
                                                depth = 10, precomputed_wavenumbers = true)

    model_base_case = NonhydrostaticModel(; grid, 
                                            boundary_conditions = (u = FieldBoundaryConditions(top = wind_stress_boundary_conditions.u),
                                                                   v = FieldBoundaryConditions(top = wind_stress_boundary_conditions.v)))

    set!(model_base_case, u = (x, y, z) -> 0.2 * (1 + randn() * 0.01), v = (x, y, z) -> randn() * 0.01)

    for n in 1:10
        time_step!(model_base_case, 0.3)
    end

    model = NonhydrostaticModel(; grid, 
                                  boundary_conditions = (u = FieldBoundaryConditions(top = wind_stress_boundary_conditions.u),
                                                         v = FieldBoundaryConditions(top = wind_stress_boundary_conditions.v)),
                                  stokes_drift)

    set!(model, u = (x, y, z) -> 0.2 * (1 + randn() * 0.01), v = (x, y, z) -> randn() * 0.01)

    for n in 1:10
        time_step!(model, 0.3)
    end

    # it instantiates and does something...
    @test !any(interior(model.velocities.u) .== interior(model_base_case.velocities.u))
end
