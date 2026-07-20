using Walrus.SurfaceFluxModel: OceanAtmosphereBoundaryConditions, 
                               PrescribedAtmosphericState, 
                               EmpiricalDownwellingLongwave,
                               AugustRocheMagnusVapourPressure,
                               heat_exchange_coefficient,
                               SimilarityTheoryInterface

@testset "Surface flux (heat and momentum)" begin
   grid = RectilinearGrid(arch; extent = (1, 1, 10), size = (1, 1, 10))
    #####
    ##### Test radiative heating/cooling (no sensible or latent flux)
    #####

    atmosphere = PrescribedAtmosphericState(wind_speed = 0.0, wind_direction = 90.0, 
                                            temperature = -272.15, 
                                            downwelling_longwave = EmpiricalDownwellingLongwave(; b = 0.0, a = 1.0, α = 0.0, γ = 0.0))

    U, V, Q = OceanAtmosphereBoundaryConditions(grid, atmosphere; 
                                                stephan_boltzman_constant = 1.0, 
                                                air_reference_density = 0.0, 
                                                water_reference_density = 1.0, 
                                                water_specific_heat_capacity = 1.0,
                                                ocean_emissivity = 0.0,
                                                controler = :T)

    model = NonhydrostaticModel(; grid, 
                                  tracers = :T,
                                  timestepper = :QuasiAdamsBashforth2,
                                  boundary_conditions = (; T = FieldBoundaryConditions(top = Q)))

    set!(model, T = -273.15)
    
    # 1 W heating -> 1 K m / s for 1s with Δz = 1m -> -272.15 K
    time_step!(model, 1)

    @test Array(interior(model.tracers.T, 1, 1, 10))[1] ≈ -272.15

    #####
    ##### Test sensible flux
    #####

    U, V, Q = OceanAtmosphereBoundaryConditions(grid, atmosphere; 
                                                stephan_boltzman_constant = 0.0, 
                                                air_reference_density = 1.0, 
                                                water_reference_density = 1.0, 
                                                water_specific_heat_capacity = 1.0,
                                                air_specific_heat_capacity = 1.0,
                                                latent_heat_vaporisation = (args...) -> 0,
                                                controler = :T)

    model = NonhydrostaticModel(; grid, 
                                  tracers = :T,
                                  timestepper = :QuasiAdamsBashforth2,
                                  boundary_conditions = (; T = FieldBoundaryConditions(top = Q)))

    set!(model, T = -273.15, u = 1)

    Ch = CUDA.@allowscalar heat_exchange_coefficient(Q.condition.func.interface_coefficients, 1, 1, grid, clock, fields(model), atmosphere)
    
    time_step!(model, 1/Ch) # 1/Cʰ at U = 1m/s

    @test Array(interior(model.tracers.T, 1, 1, 10))[1] ≈ -272.15

    #####
    ##### Test sensible heating/cooling - it is not straight forward to come up with an anlaytical result for this
    ##### so we will just check it works, then do a budget test after
    #####

    atmosphere = PrescribedAtmosphericState(wind_speed = 1.0, wind_direction = 90.0, 
                                            temperature = 0.0, 
                                            downwelling_longwave = EmpiricalDownwellingLongwave(; b = 0.0, a = 1.0, α = 0.0, γ = 0.0))

    U, V, Q = OceanAtmosphereBoundaryConditions(grid, atmosphere; 
                                                stephan_boltzman_constant=0.0, 
                                                air_reference_density = 1.0, 
                                                water_reference_density = 1.0, 
                                                water_specific_heat_capacity = 1.0,
                                                air_specific_heat_capacity = 1.0,
                                                latent_heat_vaporisation = (args...) -> 0,
                                                controler = :T)

    model = NonhydrostaticModel(; grid, 
                                  tracers = :T,
                                  timestepper = :QuasiAdamsBashforth2,
                                  boundary_conditions = (; T = FieldBoundaryConditions(top = Q)))

    # no heat exchange when temperatures equal
    set!(model, T = 0)

    for n in 1:2
        time_step!(model, 1)
    end

    @test Array(interior(model.tracers.T, 1, 1, 10))[1] ≈ 0

    # when the water is warmer, it looses heat
    set!(model, T = 1)

    for n in 1:2
        time_step!(model, 1)
    end

    @test Array(interior(model.tracers.T, 1, 1, 10))[1] < 1

    # when the water is colder, it gains heat
    set!(model, T = -1)

    for n in 1:2
        time_step!(model, 1)
    end

    @test Array(interior(model.tracers.T, 1, 1, 10))[1] > -1

    #####
    ##### Test latent heating/cooling
    #####

    atmosphere = PrescribedAtmosphericState(wind_speed = 1.0, wind_direction = 90.0, 
                                            temperature = 0.0, 
                                            downwelling_longwave = EmpiricalDownwellingLongwave(; b = 0.0, a = 1.0, α = 0.0, γ = 0.0),
                                            air_water_mixing_ratio = AugustRocheMagnusVapourPressure()(0))

    U, V, Q = OceanAtmosphereBoundaryConditions(grid, atmosphere; 
                                                stephan_boltzman_constant=0.0, 
                                                air_reference_density = 1.0, 
                                                water_reference_density = 1.0, 
                                                water_specific_heat_capacity = 1.0,
                                                air_specific_heat_capacity = 0.0,
                                                controler = :T)

    model = NonhydrostaticModel(; grid, 
                                  tracers = :T,
                                  timestepper = :QuasiAdamsBashforth2,
                                  boundary_conditions = (; T = FieldBoundaryConditions(top = Q)))

    # no heat exchange when vapour pressure equalised
    set!(model, T = 0)

    for n in 1:2
        time_step!(model, 1)
    end

    @test Array(interior(model.tracers.T, 1, 1, 10))[1] ≈ 0

    # when the water has higher saturation pressure, it looses heat
    set!(model, T = 1)
    
    for n in 1:2
        time_step!(model, 0.1)
    end

    @test Array(interior(model.tracers.T, 1, 1, 10))[1] < 1

    # when the water has lower saturation pressure, it gains heat
    set!(model, T = -1)
    
    for n in 1:2
        time_step!(model, 0.1)
    end

    @test Array(interior(model.tracers.T, 1, 1, 10))[1] > -1

    #####
    ##### Momentum transfer
    #####
    grid = RectilinearGrid(arch; size = (2, 2, 2), extent = (2, 2, 2))

    atmosphere = PrescribedAtmosphericState(wind_speed = 0.0, wind_direction = 90.0, temperature = 0.0)

    U, V, Q = OceanAtmosphereBoundaryConditions(grid, atmosphere)

    model = NonhydrostaticModel(; grid, 
                                  boundary_conditions = (u = FieldBoundaryConditions(top = U),
                                                         v = FieldBoundaryConditions(top = V)),
                                  tracers = :T)

    for n=1:10
        time_step!(model, 1)
    end

    @test all(Array(interior(model.velocities.u)) .≈ 0) & all(Array(interior(model.velocities.v)) .≈ 0) # no wind no stress

    atmosphere = PrescribedAtmosphericState(wind_speed = 1.0, wind_direction = 90.0, temperature = 0.0)

    U, V, Q = OceanAtmosphereBoundaryConditions(grid, atmosphere)

    model = NonhydrostaticModel(; grid, 
                                  boundary_conditions = (u = FieldBoundaryConditions(top = U),
                                                         v = FieldBoundaryConditions(top = V)),
                                  tracers = :T)

    time_step!(model, 1)
    
    for n=1:10000
        time_step!(model, 1)
    end

    @test -1 <= Array(interior(model.velocities.u, 1, 1, 2))[1] < 0
    @test Array(interior(model.velocities.u, 1, 1, 1))[1] ≈ 0
    @test all(Array(interior(model.velocities.v)) .≈ 0)
end
#= TODO:
@testset "Similarity theory interface coefficients" begin
    interface = SimilarityTheoryInterface()
    
end
=#
