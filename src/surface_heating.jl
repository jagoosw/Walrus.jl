module SurfaceHeatingModel

export SurfaceHeatExchange, SurfaceHeatExchangeBoundaryCondition

using Roots

using Adapt: adapt

using Oceananigans.BoundaryConditions: FluxBoundaryCondition

using Walrus: get_value, normalise_surface_function
using Walrus.WindStressModel: WindStress, 
                              LogarithmicNeutralWind

import Adapt: adapt_structure

struct SurfaceHeatExchange{WS, AT, LH, VP, FT} <: Function
                   wind_stress :: WS
               air_temperature :: AT
      latent_heat_vaporisation :: LH
               vapour_pressure :: VP
  water_specific_heat_capacity :: FT 
                 water_density :: FT
    air_specific_heat_capacity :: FT
                   air_density :: FT
        air_water_mixing_ratio :: FT
     stephan_boltzman_constant :: FT
end

adapt_structure(to, sh::SurfaceHeatExchange) = 
    SurfaceHeatExchange(adapt(to, sh.wind_stress),
                        adapt(to, sh.air_temperature),
                        adapt(to, sh.latent_heat_vaporisation),
                        adapt(to, sh.vapour_pressure),
                        sh.water_specific_heat_capacity,
                        sh.water_density,
                        sh.air_specific_heat_capacity,
                        sh.air_density,
                        sh.air_water_mixing_ratio,
                        sh.stephan_boltzman_constant)

"""
    SurfaceHeatExchange(; wind_stress,
                          air_temperature = 18, # °C
                          latent_heat_vaporisation = EmpiricalLatentHeatVaporisation(),
                          vapour_pressure = AugustRocheMagnusVapourPressure(),
                          water_specific_heat_capacity = 3991., # J / K / kg
                          water_density = 1026., # kg / m³
                          air_specific_heat_capacity = 1003.5, # J / K / kg
                          air_density = 1.204, # kg
                          air_water_mixing_ratio = 0.001, # kg / kg
                          stephan_boltzman_constant = 5.670374419e-8) # W / K⁴

Specifies surface heat exhange in the form:
```math
Q = Qᵢᵣ + Qₛ + Qₗ,
```
where ``Qᵢᵣ`` is the heat flux due to long wave (infra red) radiation, ``Qₛ`` is
the sensible heat flux, and ``Qₗ`` is the latent heat flux. Notably the short wave
radiation flux is neglegted here as it is assumed to penatrate far enough into the 
water that it is treated separately by `HomogeneousBodyHeating` (we therefore are
also assuming that the short wave penitration is suitably short that it is neglected).

The short wave term is given by the Stephan-Boltzman equation so:
```math
Qᵢᵣ = σ(T⁴ - Tₐ⁴),
```
where ``σ`` is the Stephan-Boltzman constant, ``T`` is the ocean temperature, and ``Tₐ``
is the 2m air temperature. 

The sensible and latent heat flux's are given by the bluk parameterisations described in 
[fairall2011](@citet) and are given by:
```math
Qₛ = ρₐcₚₐCₕS(T - Tₐ),
```
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
```math
F = \\frac{Q}{\\rho_oc_{po}},
```
where ``\\rho_o`` and ``c_{po}`` are the density and specific heat capacity of water.

(Note: we will retain the Oceananigans convention that negative heat flux 
 at a top boundary increases temeprature).

Keyword Arguments
=================
- `wind_stress`: wind stress model
- `air_temperature`: the air temperature in °C as a function with signature `(x, y, t)` or a constant
- `latent_heat_vaporisation`: latent heat of vaporisation in J / kg
- `vapour_pressure`: parameterisation for saturation vapour pressure in water
- `water_specific_heat_capacity`: the specific heat capacity of water in J / K / kg
- `water_density`: water density in kg / m³
- `air_specific_heat_capacit`: the specific heat capacity of air in J / K / kg
- `air_density`: air density in kg / m³
- `air_water_mixing_ratio`: water content of air in kg / kg
- `stephan_boltzman_constant`: the Stephan-Boltzman constant in W / K⁴

Example
=======

```jldoctest
julia> using Walrus

julia> using Oceananigans

julia> wind_stress = WindStress(; reference_wind_speed = 0., reference_wind_direction = 90.)
(::WindStress{Float64, Float64, LogarithmicNeutralWind{Float64, Nothing}, Float64}) (generic function with 2 methods)

julia> surface_heat_exchange = SurfaceHeatExchange(; wind_stress)
(::SurfaceHeatExchange{WindStress{Float64, Float64, LogarithmicNeutralWind{Float64, Nothing}, Float64}, Int64, Walrus.SurfaceHeatingModel.EmpiricalLatentHeatVaporisation{Float64}, Walrus.SurfaceHeatingModel.AugustRocheMagnusVapourPressure{Float64}, Float64}) (generic function with 1 method)

julia> boundary_conditions = (; T = FieldBoundaryConditions(top = FluxBoundaryCondition(surface_heat_exchange, field_dependencies = (:T, :u, :v))))
(T = Oceananigans.FieldBoundaryConditions, with boundary conditions
├── west: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── east: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── south: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── north: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── bottom: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── top: FluxBoundaryCondition: ContinuousBoundaryFunction (::SurfaceHeatExchange{WindStress{Float64, Float64, LogarithmicNeutralWind{Float64, Nothing}, Float64}, Int64, Walrus.SurfaceHeatingModel.EmpiricalLatentHeatVaporisation{Float64}, Walrus.SurfaceHeatingModel.AugustRocheMagnusVapourPressure{Float64}, Float64}) at (Nothing, Nothing, Nothing)
└── immersed: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing),)

```
"""
function SurfaceHeatExchange(; wind_stress,
                               air_temperature = 18, # °C
                               latent_heat_vaporisation = EmpiricalLatentHeatVaporisation(),
                               vapour_pressure = AugustRocheMagnusVapourPressure(),
                               water_specific_heat_capacity = 3991., # J / K / kg
                               water_density = 1026., #kg / m³
                               air_specific_heat_capacity = 1003.5, # J / K / kg
                               air_density = 1.204, # kg
                               air_water_mixing_ratio = 0.001, # kg / kgs
                               stephan_boltzman_constant = 5.670374419e-8) # W / K⁴
            
    air_temperature = normalise_surface_function(air_temperature)

    return SurfaceHeatExchange(wind_stress, air_temperature,
                               latent_heat_vaporisation, vapour_pressure,
                               water_specific_heat_capacity, water_density,
                               air_specific_heat_capacity, air_density, air_water_mixing_ratio,
                               stephan_boltzman_constant)
end

"""
    SurfaceHeatExchange(; wind_stress,
                          air_temperature = 18, # °C
                          latent_heat_vaporisation = EmpiricalLatentHeatVaporisation(),
                          vapour_pressure = AugustRocheMagnusVapourPressure(),
                          water_specific_heat_capacity = 3991., # J / K / kg
                          water_density = 1026., # kg / m³
                          air_specific_heat_capacity = 1003.5, # J / K / kg
                          air_density = 1.204, # kg
                          air_water_mixing_ratio = 0.001, # kg / kg
                          stephan_boltzman_constant = 5.670374419e-8) # W / K⁴

A convenience constructor returning `SurfaceHeatExchange` as a boundary condition

Keyword Arguments
=================
- `wind_stress`: wind stress model
- `air_temperature`: the air temperature in °C as a function with signature `(x, y, t)` or a constant
- `latent_heat_vaporisation`: latent heat of vaporisation in J / kg
- `vapour_pressure`: parameterisation for saturation vapour pressure in water
- `water_specific_heat_capacity`: the specific heat capacity of water in J / K / kg
- `water_density`: water density in kg / m³
- `air_specific_heat_capacit`: the specific heat capacity of air in J / K / kg
- `air_density`: air density in kg / m³
- `air_water_mixing_ratio`: water content of air in kg / kg
- `stephan_boltzman_constant`: the Stephan-Boltzman constant in W / K⁴

Example
=======

```jldoctest
julia> using Walrus

julia> wind_stress = WindStress(; reference_wind_speed = 0., reference_wind_direction = 90.)
(::WindStress{Float64, Float64, LogarithmicNeutralWind{Float64, Nothing}, Float64}) (generic function with 2 methods)

julia> surface_heat_exchange = SurfaceHeatExchangeBoundaryCondition(; wind_stress)
FluxBoundaryCondition: DiscreteBoundaryFunction with (::SurfaceHeatExchange{WindStress{Float64, Float64, LogarithmicNeutralWind{Float64, Nothing}, Float64}, Int64, Walrus.SurfaceHeatingModel.EmpiricalLatentHeatVaporisation{Float64}, Walrus.SurfaceHeatingModel.AugustRocheMagnusVapourPressure{Float64}, Float64})

```
"""
function SurfaceHeatExchangeBoundaryCondition(; wind_stress,
                                                air_temperature = 18, # °C
                                                latent_heat_vaporisation = EmpiricalLatentHeatVaporisation(),
                                                vapour_pressure = AugustRocheMagnusVapourPressure(),
                                                water_specific_heat_capacity = 3991., # J / K / kg
                                                water_density = 1026., #kg / m³
                                                air_specific_heat_capacity = 1003.5, # J / K / kg
                                                air_density = 1.204, # kg
                                                air_water_mixing_ratio = 0.001, # kg / kg
                                                stephan_boltzman_constant = 5.670374419e-8) # W / K⁴
            
    air_temperature = normalise_surface_function(air_temperature)

    surface_heat_exchange =  SurfaceHeatExchange(wind_stress, air_temperature,
                                                 latent_heat_vaporisation, vapour_pressure,
                                                 water_specific_heat_capacity, water_density,
                                                 air_specific_heat_capacity, air_density, air_water_mixing_ratio,
                                                 stephan_boltzman_constant)

    return FluxBoundaryCondition(surface_heat_exchange, discrete_form=true)
end

@inline function Cʰ(drag_coefficient::LogarithmicNeutralWind, wind_speed)
    κ = drag_coefficient.monin_obukhov_stability_length
    ν = drag_coefficient.air_kinematic_viscosity
    α = drag_coefficient.charnock_coefficient
    g = drag_coefficient.gravity_acceleration

    Cd = drag_coefficient(wind_speed)

    u′ = √(Cd * wind_speed)

    z₀ = α * u′^2 / g

    Rᵣ = u′ * z₀ / ν

    zₒₜ = min(1.15e-4, 5.5e-5 * Rᵣ ^ -0.6)

    result = κ  / log(10/zₒₜ) * √Cd 
  
    return ifelse(isfinite(result), result, 0)
end

# parameterisation for vapour pressure with default coefficients from [alduchov1996](@citet).
# THIS MIGHT BE WRONG!!
@kwdef struct AugustRocheMagnusVapourPressure{FT}
    a :: FT = 0.61094
    b :: FT = 17.625
    c :: FT = 243.04
end

@inline (q::AugustRocheMagnusVapourPressure)(T) = q.a * exp(q.b * T / (T + q.c)) / 1000

# parameterisation for latent heat of vaporisation for water [yu2019](@citet)
@kwdef struct EmpiricalLatentHeatVaporisation{FT}
    a :: FT = 2.501
    b :: FT = 0.00237
end

@inline (L::EmpiricalLatentHeatVaporisation)(T) = (L.a - L.b * T) * 10^6 # J / kg

#####
##### Cʰ parameterisations
#####

@inline function (interface::SurfaceHeatExchange)(i, j, grid, clock, model_fields)
    σ   = interface.stephan_boltzman_constant
    ρᵃ  = interface.air_density
    cₚᵃ = interface.air_specific_heat_capacity
    qₐ  = interface.air_water_mixing_ratio

    ρʷ  = interface.water_density
    cₚʷ = interface.water_specific_heat_capacity

    t = clock.time

    u = @inbounds model_fields.u[i, j, grid.Nz]
    v = @inbounds model_fields.v[i, j, grid.Nz]
    T = @inbounds model_fields.T[i, j, grid.Nz]

    air_temperature = get_value(interface.air_temperature, i, j, grid, clock)

    wind_speed = get_value(interface.wind_stress.reference_wind_speed, i, j, grid, clock)
    wind_direction = get_value(interface.wind_stress.reference_wind_direction, i, j, grid, clock)

    uʷ = - wind_speed * sind(wind_direction)
    vʷ = - wind_speed * cosd(wind_direction)

    relative_speed = √((uʷ - u)^2 + (vʷ - v)^2)

    cʰ = Cʰ(interface.wind_stress.drag_coefficient, relative_speed)
    
    radiative_cooling = σ * ((273.15 + T) ^ 4 - (273.15 + air_temperature) ^ 4) # J / m² / s

    sensible_cooling = relative_speed * ρᵃ * cₚᵃ * cʰ * (T - air_temperature) # (m / s) (kg / m³) (J / kg / K) (1) (k) -> (1 / s) (1 / m²) (J) -> J / m² / s

    q = interface.vapour_pressure(T)
    L = interface.latent_heat_vaporisation(T)

    latent_cooling = relative_speed * ρᵃ * L * cʰ * (q - qₐ) # (m / s) (kg / m³) (J / kg) (1) (kg / kg) -> (m / s) (m³) (J) -> J / m² / s

    return (radiative_cooling + sensible_cooling + latent_cooling) / (ρʷ * cₚʷ) # (J / m² / s) / ((kg / m³) (J / kg / K)) -> (J / m² / s) / ( J / m³ / K)) -> K m / s
end

end # module
