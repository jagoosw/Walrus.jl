using Oceananigans.BoundaryConditions: FluxBoundaryCondition, BoundaryCondition, DiscreteBoundaryFunction

struct OceanAtmosphereBoundary{IC, AS, WD, AD, WC, AC, VP, LH, FT} <: Function
          interface_coefficients :: IC
                atmosphere_state :: AS

         water_reference_density :: WD
           air_reference_density :: AD
    water_specific_heat_capacity :: WC
      air_specific_heat_capacity :: AC
                 vapour_pressure :: VP
        latent_heat_vaporisation :: LH
       stephan_boltzman_constant :: FT
                ocean_emissivity :: FT
end

Adapt.adapt_structure(to, boundary::OceanAtmosphereBoundary) =
    OceanAtmosphereBoundary(adapt(to, boundary.interface_coefficients),
                            adapt(to, boundary.atmosphere_state),
                            adapt(to, boundary.water_reference_density),
                            adapt(to, boundary.air_reference_density),
                            adapt(to, boundary.water_specific_heat_capacity),
                            adapt(to, boundary.air_specific_heat_capacity),
                            adapt(to, boundary.vapour_pressure),
                            adapt(to, boundary.latent_heat_vaporisation),
                            boundary.stephan_boltzman_constant,
                            boundary.ocean_emissivity)

function OceanAtmosphereBoundaryConditions(grid, atmosphere_state; 
                                           interface_coefficients = SimilarityTheoryInterface(grid),
                                           water_reference_density = 1026.0, # TODO: make this a function of temperature and salinity
                                           air_reference_density = 1.225,
                                           water_specific_heat_capacity = 3991., # J / K / kg,  TODO: make this a function of temperature and salinity
                                           air_specific_heat_capacity = 1003.5, # J / K / kg
                                           vapour_pressure = AugustRocheMagnusVapourPressure(),
                                           latent_heat_vaporisation = EmpiricalLatentHeatVaporisation(),
                                           stephan_boltzman_constant = 5.670374419e-8, # W / K⁴
                                           ocean_emissivity = 0.97)

    boundary = OceanAtmosphereBoundary(interface_coefficients, atmosphere_state,
                                       water_reference_density, air_reference_density, 
                                       water_specific_heat_capacity, air_specific_heat_capacity, 
                                       vapour_pressure, latent_heat_vaporisation,
                                       stephan_boltzman_constant, ocean_emissivity)

    u = FluxBoundaryCondition(boundary; parameters = Val(:u), discrete_form=true)
    v = FluxBoundaryCondition(boundary; parameters = Val(:v), discrete_form=true)
    T = FluxBoundaryCondition(boundary; parameters = Val(:T), discrete_form=true)
    # TODO: add evaporation and proper calculations for other scalars like CO₂

    return (; u, v, T)
end

@inline function (boundary::OceanAtmosphereBoundary)(i, j, grid, clock, model_fields, ::Val{:u})
    ρₐ = boundary.air_reference_density
    ρₒ = boundary.water_reference_density

    uʷ = relative_x_wind(boundary.atmosphere_state, i, j, grid, clock, model_fields)
    U  = relative_wind_speed(boundary.atmosphere_state, i, j, grid, clock, model_fields)

    Cd = drag_coefficient(boundary.interface_coefficients, i, j, grid, clock, model_fields, boundary.atmosphere_state)

    stress_velocity = ρₐ / ρₒ * Cd * U

    return - stress_velocity * uʷ
end

@inline function (boundary::OceanAtmosphereBoundary)(i, j, grid, clock, model_fields, ::Val{:v})
    ρₐ = boundary.air_reference_density
    ρₒ = boundary.water_reference_density

    vʷ = relative_y_wind(boundary.atmosphere_state, i, j, grid, clock, model_fields)
    U  = relative_wind_speed(boundary.atmosphere_state, i, j, grid, clock, model_fields)

    Cd = drag_coefficient(boundary.interface_coefficients, i, j, grid, clock, model_fields, boundary.atmosphere_state)

    stress_velocity = ρₐ / ρₒ * Cd * U

    return - stress_velocity * vʷ
end

@inline function (boundary::OceanAtmosphereBoundary)(i, j, grid, clock, model_fields, ::Val{:T})
    ρₐ  = boundary.air_reference_density
    ρₒ  = boundary.water_reference_density
    σ   = boundary.stephan_boltzman_constant
    cₚᵃ = boundary.air_specific_heat_capacity
    cₚʷ = boundary.water_specific_heat_capacity
    ϵ   = boundary.ocean_emissivity
    
    T = @inbounds model_fields.T[i, j, grid.Nz]

    qₐ = air_water_mixing_ratio(boundary.atmosphere_state, i, j, grid, clock, model_fields)
    θ  = temperature(boundary.atmosphere_state, i, j, grid, clock, model_fields)
    U  = relative_wind_speed(boundary.atmosphere_state, i, j, grid, clock, model_fields)

    q = boundary.vapour_pressure(T)
    L = boundary.latent_heat_vaporisation(T)

    Ch = heat_exchange_coefficient(boundary.interface_coefficients, i, j, grid, clock, model_fields, boundary.atmosphere_state)

    radiative_cooling = ϵ * σ * (273.15 + T) ^ 4 - downwelling_longwave(boundary.atmosphere_state, i, j, grid, clock, model_fields, boundary) # J / m² / s

    sensible_cooling = ρₐ * cₚᵃ * Ch * (T - θ) * U # (m / s) (kg / m³) (J / kg / K) (1) (k) -> (1 / s) (1 / m²) (J) -> J / m² / s

    latent_cooling = ρₐ * L * Ch * (q - qₐ) * U # (m / s) (kg / m³) (J / kg) (1) (kg / kg) -> (m / s) (m³) (J) -> J / m² / s

    return (radiative_cooling + sensible_cooling + latent_cooling) / (ρₒ * cₚʷ) # (J / m² / s) / ((kg / m³) (J / kg / K)) -> (J / m² / s) / ( J / m³ / K)) -> K m / s
end


#####
##### parameterisation for vapour pressure with default coefficients from [alduchov1996](@citet).
#####

@kwdef struct AugustRocheMagnusVapourPressure{FT}
   e0 :: FT = 6.122
    p :: FT = 1013.0
    a :: FT = 0.61094
    b :: FT = 17.625
    c :: FT = 243.04
end

@inline function (q::AugustRocheMagnusVapourPressure)(T)
  es = water_vapour_pressure(q, T)

  return q.a * es / (q.p - (1 - q.a) * es)
end

@inline water_vapour_pressure(q::AugustRocheMagnusVapourPressure, T) = 
    q.e0 * exp(q.b * T / (T + q.c))


#####
##### parameterisation for latent heat of vaporisation for water [yu2019](@citet)
#####

@kwdef struct EmpiricalLatentHeatVaporisation{FT}
    a :: FT = 2.501
    b :: FT = 0.00237
end

@inline (L::EmpiricalLatentHeatVaporisation)(T) = (L.a - L.b * T) * 10^6 # J / kg


#####
##### update coefficients
#####

function update_boundary_condition!(bc::BoundaryCondition{<:Any, <:DiscreteBoundaryFunction{<:Any, <:OceanAtmosphereBoundary}}, ::Val{:top}, field, model)
    interface = bc.condition.func.interface_coefficients
    atmosphere = bc.condition.func.atmosphere_state

    update_interface!(interface, model, atmosphere)

    return nothing
end