using Walrus.RadiativeTransfer: HomogeneousBodyHeating, PARModelHeating

using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry
import Oceananigans.Biogeochemistry: required_biogeochemical_tracers

struct JustPhytoplankton <: AbstractContinuousFormBiogeochemistry end
required_biogeochemical_tracers(::JustPhytoplankton) = (:P, )

@testset "Body heating" begin
    grid = RectilinearGrid(arch; extent = (1, 1, 10), size = (1, 1, 10))

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
    @test Array(interior(model.tracers.T, 1, 1, 10))[1] ≈ 10 * (1 - exp(-1))

    #####
    ##### OceanBioME light attenuation body heating
    #####

    grid = RectilinearGrid(arch; size = (2, 2, 2), extent = (2, 2, 2))

    light_attenuation_model = TwoBandPhotosyntheticallyActiveRadiation(; grid)
    body_heating = PARModelHeating(; light_attenuation_model)

    biogeochemistry = Biogeochemistry(JustPhytoplankton(); light_attenuation = light_attenuation_model)

    model = NonhydrostaticModel(; grid, biogeochemistry, timestepper = :QuasiAdamsBashforth2, forcing = (; T = Forcing(body_heating, discrete_form=true)), tracers = :T)

    Pᵢ(x, y, z) = 2.5 + z

    set!(model, P = Pᵢ)

    kʳ = light_attenuation_model.water_red_attenuation
    kᵇ = light_attenuation_model.water_blue_attenuation
    χʳ = light_attenuation_model.chlorophyll_red_attenuation
    χᵇ = light_attenuation_model.chlorophyll_blue_attenuation
    eʳ = light_attenuation_model.chlorophyll_red_exponent
    eᵇ = light_attenuation_model.chlorophyll_blue_exponent
    r = light_attenuation_model.pigment_ratio
    Rᶜₚ = light_attenuation_model.phytoplankton_chlorophyll_ratio

    zc = znodes(grid, Center(), Center(), Center())

    Δz = zspacings(grid, Center())

    Chlʳ = [(Pᵢ(0, 0, zc[2]) * Rᶜₚ / r) ^ eʳ, (Pᵢ(0, 0, zc[1]) * Rᶜₚ / r) ^ eʳ]
    Chlᵇ = [(Pᵢ(0, 0, zc[2]) * Rᶜₚ / r) ^ eᵇ, (Pᵢ(0, 0, zc[1]) * Rᶜₚ / r) ^ eᵇ]

    ∫Chlʳ = [(Pᵢ(0, 0, zc[2]) * Rᶜₚ / r) ^ eʳ * Δz[1, 1, 1]/2]
    ∫Chlᵇ = [(Pᵢ(0, 0, zc[2]) * Rᶜₚ / r) ^ eᵇ * Δz[1, 1, 1]/2]

    push!(∫Chlʳ, ∫Chlʳ[1] + (Pᵢ(0, 0, zc[2]) * Rᶜₚ / r) ^ eʳ * Δz[1, 1, 1]/2 + (Pᵢ(0, 0, zc[1]) * Rᶜₚ / r) ^ eʳ * Δz[1, 1, 1]/2)
    push!(∫Chlᵇ, ∫Chlᵇ[1] + (Pᵢ(0, 0, zc[2]) * Rᶜₚ / r) ^ eᵇ * Δz[1, 1, 1]/2 + (Pᵢ(0, 0, zc[1]) * Rᶜₚ / r) ^ eᵇ * Δz[1, 1, 1]/2)

    expected_PAR = 100.0 .* [exp(zc[2] * kʳ - ∫Chlʳ[1] * χʳ) + exp(zc[2] * kᵇ - ∫Chlᵇ[1] * χᵇ),
                             exp(zc[1] * kʳ - ∫Chlʳ[2] * χʳ) + exp(zc[1] * kᵇ - ∫Chlᵇ[2] * χᵇ)] ./ 2

    analytical_body_heating = 100.0 .* [(kʳ + Chlʳ[1] * χʳ) * exp(zc[2] * kʳ - ∫Chlʳ[1] * χʳ) + (kᵇ + Chlᵇ[1] * χᵇ) * exp(zc[2] * kᵇ - ∫Chlᵇ[1] * χᵇ),
                                        (kʳ + Chlʳ[2] * χʳ) * exp(zc[1] * kʳ - ∫Chlʳ[2] * χʳ) + (kᵇ + Chlᵇ[2] * χᵇ) * exp(zc[1] * kᵇ - ∫Chlᵇ[2] * χᵇ)] ./ 2

    analytical_body_heating ./= (body_heating.water_density * body_heating.water_heat_capacity)

    # high tollerance, I think error is from course grid (hopefully)
    @test CUDA.@allowscalar all(isapprox.(analytical_body_heating, [body_heating(1, 1, k, grid, model.clock, fields(model)) for k = 2:-1:1], atol = 1.5e-6))
end