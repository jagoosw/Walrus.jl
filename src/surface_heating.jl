module SurfaceHeatingModel

export AirSeaInterface

using Walrus: ReturnValue, display_input
using Walrus.WindStressModel: WindStress, LogarithmicNeutralWind

struct SurfaceHeatExchange{CD, AT, LH, VP, FT} <: Function
              drag_coefficient :: CD
               air_temperature :: AT
      latent_heat_vaporisation :: LH
               vapour_pressure :: VP
  water_specific_heat_capacity :: FT 
                 water_density :: FT
    air_specific_heat_capacity :: FT
                   air_density :: FT
     stephan_boltzman_constant :: FT
end

"""
    SurfaceHeatExchange(; drag_coefficient = LogarithmicNeutralWind(),
                          air_temperature = 18, # °C
                          latent_heat_vaporisation = EmpiricalLatentHeatVaporisation(),
                          vapour_pressure = AugustRocheMagnusVapourPressure(),
                          water_specific_heat_capacity = 3991., # J / K / kg
                          water_density = 1026., # kg / m³
                          air_specific_heat_capacity = 1003.5, # J / K / kg
                          air_density = 1.204, # kg
                          stephan_boltzman_constant = 5.670374419e-8) # W / K⁴

Specifies surface heat exhange in the form:
``Q = Qᵢᵣ + Qₛ + Qₗ``, 
where ``Qᵢᵣ`` is the heat flux due to long wave (infra red) radiation, ``Qₛ`` is
the sensible heat flux, and ``Qₗ`` is the latent heat flux. Notably the short wave
radiation flux is neglegted here as it is assumed to penatrate far enough into the 
water that it is treated separately by ``HomogeneousBodyHeating`` (we therefore are
also assuming that the short wave penitration is suitably short that it is neglected).

### Q is J / m² / s

The short wave term is given by the Stephan-Boltzman equation so:
``Qᵢᵣ = σ(T⁴ - Tₐ⁴)``,
where ``σ`` is the Stephan-Boltzman constant, ``T`` is the ocean temperature, and ``Tₐ``
is the 2m air temperature. 

The sensible and latent heat flux's are given by the bluk parameterisations described in 
[fairall2011](@citet) and are given by:
``Qₛ = ρₐcₚₐCₕS(T - Tₐ)``, 
and ``Qₗ = ρₐLₑCₕS(q(T) - qₐ)``, 
where ``ρₐ`` is the density of the air, ``cₚₐ`` is the specific heat capacity of air, 
``Cₕ`` is the heat transfer coefficient, ``S`` is the wind speed, ``Lₑ`` is the latent
heat of vaporizaion parameterised by `latent_heat_vaporisation`, ``q(T)`` is the 
saturation vapour pressure of water and is parameterised by `vapour_pressure`.

The heat transfer coefficients are given by the same parameterisation as the `WindStress`.
For example for the `LogarithmicNeutralWind` the transfer coefficient is given as:
``C_h = \\frac{\\kappa}{\\log{\\frac{2}{z_0}}}\\frac{\\kappa}{\\log{\\frac{2}{z_{ot}}}}``,
where ``z_{ot}`` is the sclar roughness parameter given by:
``z_{ot} = \\min\\left(1.15\\cdot10^{-4}, 5.5\\cdot10^{-5}R_r^{-0.6}\\right)``,
where ``R_r`` is the roughness reynolds number given as ``R_r = \\frac{u\\star z_0}{\\nu}``
where ``\\nu`` is the kinematic viscosity of air.

The heat flux is then given by:
``F = \\frac{Q}{\\rho_oc_{po}}``,
where ``\\rho_o`` and ``c_{po}`` are the density and specific heat capacity of water.

(Note: we will retain the Oceananigans convention that negative heat flux 
 at a top boundary increases temeprature).
"""
function SurfaceHeatExchange(; drag_coefficient = LogarithmicNeutralWind(),
                               air_temperature = 18, # °C
                               latent_heat_vaporisation = EmpiricalLatentHeatVaporisation(),
                               vapour_pressure = AugustRocheMagnusVapourPressure(),
                               water_specific_heat_capacity = 3991., # J / K / kg
                               water_density = 1026., #kg / m³
                               air_specific_heat_capacity = 1003.5, # J / K / kg
                               air_density = 1.204, # kg
                               stephan_boltzman_constant = 5.670374419e-8) # W / K⁴
            
    isa(air_temperature, Function) || (air_temperature = ReturnValue(air_temperature))

    return SurfaceHeatExchange(drag_coefficient, air_temperature,
                               latent_heat_vaporisation, vapour_pressure,
                               water_specific_heat_capacity, water_density,
                               air_specific_heat_capacity, air_density,
                               stephan_boltzman_constant)
end

@inline function Cʰ(drag_coefficient::LogarithmicNeutralWind, wind_speed)
    κ = drag_coefficient.monin_obukhov_stability_length
    ν = drag_coefficient.air_kinematic_viscosity
    aᶜ = drag_coefficient.charnock_coefficient
    b = drag_coefficient.gravity_wave_coefficient
    g = drag_coefficient.gravity_acceleration

    params = (; κ, ν, aᶜ, b, g, wind_speed)

    ū = find_zero(velocity_roughness_length_roots, 1, p = params)
    z₀ = velocity_roughness_length(ū, params)

    wind_speed == 0 && (z₀ = Inf)

    Rᵣ = ū * z₀ / params.ν
    zₒₜ = min(1.15e-4, 5.5e-5 * Rᵣ ^ -0.6)

    return params.κ ^ 2 / (log(2/z₀) * log(2/zₒₜ))
end

# parameterisation for vapour pressure with default coefficients from [alduchov1996](@citet).
@kwdef struct AugustRocheMagnusVapourPressure{FT}
    a :: FT = 0.61094
    b :: FT = 17.625
    c :: FT = 243.04
end

@inline (q::AugustRocheMagnusVapourPressure)(T) = q.a * exp(q.b * T / (T + q.c)) / 1000 # 

# parameterisation for latent heat of vaporisation for water [yu2019](@citet)
@kwdef struct EmpiricalLatentHeatVaporisation{FT}
    a :: FT = 2.501
    b :: FT = 0.00237
end

@inline (L::EmpiricalLatentHeatVaporisation)(T) = (L.a - L.b * T) * 10^6 # J / kg

#####
##### Cʰ parameterisations
#####

@inline function (interface::SurfaceHeatExchange)(x, y, t, T, u, v)
    σ   = interface.stephan_boltzman_constant
    ρᵃ  = interface.air_density
    cₚᵃ = interface.air_specific_heat_capacity
    qₐ  = interface.air_water_mixing_ratio

    ρʷ  = interface.water_density
    cₚʷ = interface.water_specific_heat_capacity

    air_temperature = interface.air_temperature(x, y, t)

    wind_speed = interface.wind_stress.reference_wind_speed(x, y, t)

    relative_speed = wind_speed - sqrt(u^2 + v^2)

    cʰ = Cʰ(interface.wind_stress.drag_coefficient, relative_speed)
    
    radiative_cooling = σ * ((273.15 + T) ^ 4 - (273.15 + air_temperature) ^ 4)

    sensible_cooling = relative_speed * ρᵃ * cₚᵃ * cʰ * (T - air_temperature)

    q = interface.vapour_mixing_ration(T)
    L = interface.latent_heat_vaporisation(T)

    latent_cooling = relative_speed * ρᵃ * cₚᵃ * L * (q - qₐ)

    return (radiative_cooling + sensible_cooling + latent_cooling) / (ρʷ * cₚʷ)
end

end # module