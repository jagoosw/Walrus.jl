@testset "Surface heat exchange" begin
    grid = RectilinearGrid(extent = (1, 1, 10), size = (1, 1, 10))

    #####
    ##### Test body heating
    #####

    body_heating = HomogeneousBodyHeating(; surface_flux = (x, y, t) -> 1,
                                            water_attenuation_coefficient = 1.,
                                            water_density = 1.,
                                            water_heat_capacity = 1.)

    model = NonhydrostaticModel(; grid, forcing = (; T = Forcing(body_heating, discrete_form=true)), tracers = :T)

    time_step!(model, 10)

    # 1 W/m² for 10s -> 10J, top meter absorbs 10J * (exp(0) - exp(- 1)) ≈ 6.3212 J
    # Density and specific head are 1 so the water should have increased in temp by ~6.3212 K
    @test model.tracers.T[1, 1, 10] ≈ 10 * (1 - exp(-1))

    # TODO: change this test so it is actually correct when spacing is not 1 in all dimensions

    #####
    ##### Test radiative heating/cooling (no sensible or latent flux)
    #####

    wind_stress_boundary_conditions = WindStressBoundaryConditions(; reference_wind_speed = 0., reference_wind_direction = 90.)

    surface_heat_exchange = SurfaceHeatExchangeBoundaryCondition(; wind_stress = wind_stress_boundary_conditions.u.condition.func,
                                                                   air_temperature = -272.15,
                                                                   air_density = 0.,
                                                                   stephan_boltzman_constant = 1.,
                                                                   water_density = 1.,
                                                                   water_specific_heat_capacity = 1.)

    model = NonhydrostaticModel(; grid, 
                                  tracers = :T,
                                  boundary_conditions = (; T = FieldBoundaryConditions(top = surface_heat_exchange)))

    set!(model, T = -273.15)
    
    # 1 W heating -> 1 K m / s for 1s with Δz = 1m -> -272.15 K
    time_step!(model, 1)

    @test model.tracers.T[1, 1, 10] ≈ -272.15

    #####
    ##### Test sensible flux
    #####

    wind_stress_boundary_conditions = WindStressBoundaryConditions(; reference_wind_speed = 0., reference_wind_direction = 90.)

    surface_heat_exchange = SurfaceHeatExchangeBoundaryCondition(; wind_stress = wind_stress_boundary_conditions.u.condition.func,
                                                                   air_temperature = -272.15,
                                                                   air_density = 0.,
                                                                   stephan_boltzman_constant = 1.,
                                                                   water_density = 1.,
                                                                   water_specific_heat_capacity = 1.)

    model = NonhydrostaticModel(; grid, 
                                  tracers = :T,
                                  boundary_conditions = (; T = FieldBoundaryConditions(top = surface_heat_exchange)))

    set!(model, T = -273.15)
    
    # 1 W heating -> 1 K m / s for 1s with Δz = 1m -> -272.15 K
    time_step!(model, 1)

    @test model.tracers.T[1, 1, 10] ≈ -272.15

    #####
    ##### Test sensible heating/cooling - it is not straight forward to come up with an anlaytical result for this
    ##### so we will just check it works, then do a budget test after
    #####

    wind_stress_boundary_conditions = WindStressBoundaryConditions(; reference_wind_speed = 1., reference_wind_direction = 90.)

    surface_heat_exchange = SurfaceHeatExchangeBoundaryCondition(; wind_stress = wind_stress_boundary_conditions.u.condition.func,
                                                                   air_temperature = 0,
                                                                   air_density = 1.,
                                                                   air_specific_heat_capacity = 1.,
                                                                   stephan_boltzman_constant = 0.,
                                                                   water_density = 1.,
                                                                   water_specific_heat_capacity = 1.,
                                                                   latent_heat_vaporisation = (T) -> 0)

    model = NonhydrostaticModel(; grid, 
                                  tracers = :T,
                                  boundary_conditions = (; T = FieldBoundaryConditions(top = surface_heat_exchange)))

    # no heat exchange when temperatures equal
    set!(model, T = 0)
    
    time_step!(model, 1)

    @test model.tracers.T[1, 1, 10] ≈ 0

    # when the water is warmer, it looses heat
    set!(model, T = 1)
    
    time_step!(model, 1)

    @test model.tracers.T[1, 1, 10] < 1

    # when the water is colder, it gains heat
    set!(model, T = -1)
    
    time_step!(model, 1)

    @test model.tracers.T[1, 1, 10] > -1

    #####
    ##### Test latent heating/cooling
    #####

    wind_stress_boundary_conditions = WindStressBoundaryConditions(; reference_wind_speed = 1., reference_wind_direction = 90.)

    surface_heat_exchange = SurfaceHeatExchangeBoundaryCondition(; wind_stress = wind_stress_boundary_conditions.u.condition.func,
                                                                   air_temperature = 0,
                                                                   air_density = 1.,
                                                                   air_specific_heat_capacity = 0.,
                                                                   air_water_mixing_ratio = 0.00061094,
                                                                   stephan_boltzman_constant = 0.,
                                                                   water_density = 1.,
                                                                   water_specific_heat_capacity = 1.)

    model = NonhydrostaticModel(; grid, 
                                  tracers = :T,
                                  boundary_conditions = (; T = FieldBoundaryConditions(top = surface_heat_exchange)))

    # no heat exchange when vapour pressure equalised
    set!(model, T = 0)
    
    time_step!(model, 1)

    @test model.tracers.T[1, 1, 10] ≈ 0

    # when the water has higher saturation pressure, it looses heat
    set!(model, T = 1)
    
    time_step!(model, 1)

    @test model.tracers.T[1, 1, 10] < 1

    # when the water has lower saturation pressure, it gains heat
    set!(model, T = -1)
    
    time_step!(model, 1)

    @test model.tracers.T[1, 1, 10] > -1
end
