module WindStressModel

using Walrus: ReturnValue, display_input

export WindStress, WindStressBoundaryConditions, LogarithmicNeutralWind

using Oceananigans.BuoyancyModels: g_Earth

import Base: summary, show

struct WindStress{WS, WD, DC, FT} <: Function
      reference_wind_speed :: WS
  reference_wind_direction :: WD
          drag_coefficient :: DC
               air_density :: FT
             water_density :: FT
end

"""
    WindStress(; reference_wind_speed, 
                 reference_wind_direction,
                 drag_coefficient = LogarithmicNeutralWind(), 
                 air_density = 1.225, 
                 water_density = 1026.)

Returns a wind stress model where the stress is given by,
``\\frac{\\tau}{\\rho_o} = \\rho_aC_dSU_{x/y}``,
where ``\\rho_o`` is the water density, ``\\rho_a`` is the air density,
``C_d`` is the drag coefficient, ``U_{x/y}`` are the x and y components of 
relative wind speed, and ``S=\\sqrt{U_x^2+U_y^2}``.

``C_d`` is calculated from a parameterisation, by default this is a "log neutral" wind
parameterisation with velocity roughness length parameterisaion like [smith1988](@citet).

In the default configuration this is the same as described in [fairall2011](@citet).

Keyword Arguments
=================

- `reference_wind_speed` (required): a function returning the (10m neutral) wind speed in the form `reference_wind_speed(x, y, t)` or single value
- `reference_wind_direction` (required): a function returning the (10m neutral) wind direction in the form `reference_wind_direction(x, y, t)` or single value
- `drag_coefficient`: the drag coefficient parameterisation
- `air_density`: air density in kg/m³ 
- `water_density`: water density in kg/m³ 

Example
=======

```jldoctest
julia> using Walrus: WindStress

julia> using Oceananigans

julia> reference_wind_speed = 0.1
0.1

julia> reference_wind_direction = 0.
0.0

julia> wind_stress = WindStress(; reference_wind_speed, reference_wind_direction)
(::WindStress{Walrus.ReturnValue{Float64}, Walrus.ReturnValue{Float64}, LogarithmicNeutralWind{Float64}, Float64}) (generic function with 2 methods)

julia> boundary_conditions = (u = FieldBoundaryConditions(top = FluxBoundaryCondition(wind_stress, parameters = Val(:x))),
                              v = FieldBoundaryConditions(top = FluxBoundaryCondition(wind_stress, parameters = Val(:y))))
(u = Oceananigans.FieldBoundaryConditions, with boundary conditions
├── west: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── east: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── south: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── north: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── bottom: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── top: FluxBoundaryCondition: ContinuousBoundaryFunction (::WindStress{Walrus.ReturnValue{Float64}, Walrus.ReturnValue{Float64}, LogarithmicNeutralWind{Float64}, Float64}) at (Nothing, Nothing, Nothing)
└── immersed: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing), v = Oceananigans.FieldBoundaryConditions, with boundary conditions
├── west: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── east: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── south: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── north: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── bottom: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── top: FluxBoundaryCondition: ContinuousBoundaryFunction (::WindStress{Walrus.ReturnValue{Float64}, Walrus.ReturnValue{Float64}, LogarithmicNeutralWind{Float64}, Float64}) at (Nothing, Nothing, Nothing)
└── immersed: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing))

```
"""
function WindStress(; reference_wind_speed, 
                      reference_wind_direction,
                      drag_coefficient = LogarithmicNeutralWind(), 
                      air_density = 1.225, 
                      water_density = 1026.)
        
    isa(reference_wind_speed, Function) || (reference_wind_speed = ReturnValue(reference_wind_speed))

    isa(reference_wind_direction, Function) || (reference_wind_direction = ReturnValue(reference_wind_direction))

    return WindStress(reference_wind_speed, reference_wind_direction,
                      drag_coefficient, air_density, water_density)
end

"""
    WindStressBoundaryConditions(; reference_wind_speed, 
                                   reference_wind_direction,
                                   drag_coefficient = LogarithmicNeutralWind(), 
                                   air_density = 1.225, 
                                   water_density = 1026.)


Convenience constructor to setup `WindStress` boundary conditions.

Keyword Arguments
=================

- `reference_wind_speed` (required): a function returning the (10m neutral) wind speed in the form `reference_wind_speed(x, y, t)` or single value
- `reference_wind_direction` (required): a function returning the (10m neutral) wind direction in the form `reference_wind_direction(x, y, t)` or single value
- `drag_coefficient`: the drag coefficient parameterisation
- `air_density`: air density in kg/m³ 
- `water_density`: water density in kg/m³ 

Example
=======

```jldoctest
julia> using Walrus: WindStressBoundaryConditions

julia> wind_stress_boundary_conditions = WallStressBoundaryConditions()
(u = FluxBoundaryCondition: DiscreteBoundaryFunction (::WallStress{Float64}) with parameters Val{:x}, v = FluxBoundaryCondition: DiscreteBoundaryFunction (::WallStress{Float64}) with parameters Val{:y})

julia> boundary_conditions = (u = FieldBoundaryConditions(top = wind_stress_boundary_conditions.u),
                              v = FieldBoundaryConditions(top = wind_stress_boundary_conditions.v))
                              (u = Oceananigans.FieldBoundaryConditions, with boundary conditions
                              ├── west: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
                              ├── east: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
                              ├── south: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
                              ├── north: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
                              ├── bottom: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
                              ├── top: FluxBoundaryCondition: DiscreteBoundaryFunction (::WallStress{Float64}) with parameters Val{:x}
                              └── immersed: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing), v = Oceananigans.FieldBoundaryConditions, with boundary conditions
                              ├── west: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
                              ├── east: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
                              ├── south: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
                              ├── north: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
                              ├── bottom: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
                              ├── top: FluxBoundaryCondition: DiscreteBoundaryFunction (::WallStress{Float64}) with parameters Val{:y}
                              └── immersed: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing))

```
"""
function WindStressBoundaryConditions(; reference_wind_speed, 
                                        reference_wind_direction,
                                        drag_coefficient = LogarithmicNeutralWind(), 
                                        air_density = 1.225, 
                                        water_density = 1026.)

    wind_stress = WindStress(; reference_wind_speed, reference_wind_direction,
                               drag_coefficient, air_density, water_density)

    u = FluxBoundaryCondition(wind_stress, parameters = Val(:x))

    v = FluxBoundaryCondition(wind_stress, parameters = Val(:y))

    return (; u, v)
end

@inline function (wind_stress::WindStress)(x, y, t, u, v, ::Val{:x})
    ρₐ = wind_stress.air_density
    ρₒ = wind_stress.water_density

    wind_speed = wind_stress.reference_wind_speed(x, y, t)
    wind_direction = wind_stress.reference_wind_direction(x, y, t)

    relative_speed = wind_speed - sqrt(u^2 + v^2)

    stress_velocity = ρₐ / ρₒ * wind_speed.drag_coefficient(relative_speed) * relative_speed

    uʷ = - wind_speed * sind(wind_direction)

    return - stress_velocity * (uʷ - u)
end

@inline function (wind_stress::WindStress)(x, y, t, u, v, ::Val{:y})
    ρₐ = wind_stress.air_density
    ρₒ = wind_stress.water_density

    wind_speed = wind_stress.reference_wind_speed(x, y, t)
    wind_direction = wind_stress.reference_wind_direction(x, y, t)

    relative_speed = wind_speed - sqrt(u^2 + v^2)

    stress_velocity = ρₐ / ρₒ * wind_speed.drag_coefficient(relative_speed) * relative_speed

    vʷ = - wind_speed * cosd(wind_direction)

    return - stress_velocity * (vʷ - v)
end

summary(::WindStress) = string("Wind stress model")
show(io::IO, wind::WindStress) = println(io, summary(wind), " with:\n",
                                     " Wind speed: ", display_input(wind.reference_wind_speed), "\n",
                                     " Wind direction: ", display_input(wind.reference_wind_direction), "\n",
                                     " Drag coefficient: ", summary(wind.drag_coefficient), "\n",
                                     " Air density: ", wind.air_density, " kg/m³\n",
                                     " Water density: ", wind.water_density, " kg/m³")

"""
    LogarithmicNeutralWind(; monin_obukhov_stability_length = 0.4
                             charnock_coefficient = 0.014
                             air_kinematic_viscosity = 1.488e-5
                             gravity_wave_coefficient = 0.11
                             gravity = g_Earth)

Returns a `LogarithmicNeutralWind` parameterisation for the surface drag coefficient

``C_d`` is parameterised as,
``C_d = \\frac{\\kappa^2}{\\log{\\frac{10}{z_0}}}``,
where ``\\kappa`` is the Monin‐Obukhov stability length and ``z_0`` is the velocity 
roughness length. This is the roughness length scale which logarithmically brings 
the relative velocity to zero at the surface, i.e.
``U=\\frac{u\\star}{\\kappa}\\log\\frac{z/z_0}``,
where ``u\\star`` is the friction velocity. Additionally ``z_0`` is given as,
``z_0=b\\frac{\\nu}{u\\star} + \\frac{a_c}{g}u\\star^2``,
where ``\\nu`` is the kinematic viscosity of air and g is the acceleration of gravity.

This model itterativly solves these equations to find ``z_0``.

This parameterisaion is described in [smith1988](@citet)
"""
@kwdef struct LogarithmicNeutralWind{FT} # can't think of a good name for this
  monin_obukhov_stability_length :: FT = 0.4
            charnock_coefficient :: FT = 0.014
         air_kinematic_viscosity :: FT = 1.488e-5
        gravity_wave_coefficient :: FT = 0.11
                         gravity :: FT = g_Earth
end

@inline velocity_roughness_length_roots(ū, params) = velocity_roughness_length(ū, params) - 10 * exp(- params.wind_speed / (params.κ * ū))

@inline velocity_roughness_length(ū, params) = 0.11 * params.ν / (ū + eps(0.0)) + params.aᶜ * ū^2 / params.gravity

@inline function (dc::LogarithmicNeutralWind)(wind_speed)
    κ = dc.monin_obukhov_stability_length
    ν = dc.air_kinematic_viscosity
    aᶜ = dc.charnock_coefficient
    b = dc.gravity_wave_coefficient

    ū = find_zero(velocity_roughness_length_roots, 1, p = (; κ, ν, aᶜ, b, wind_speed))
    z₀ = velocity_roughness_length(ū, params)

    return (params.κ / log(10/z₀) ^ 2)
end

summary(::LogarithmicNeutralWind) = string("Log neutral wind drag coefficient model")
show(io::IO, cd::LogarithmicNeutralWind) = println(io, summary(cd), " with:\n",
                                                  " κ: ",  cd.monin_obukhov_stability_length, "\n",
                                                  " ν: ",  cd.air_kinematic_viscosity, " m²/s\n",
                                                  " aᶜ: ", cd.charnock_coefficient, "\n",
                                                  " g: " , cd.gravity, "m/s²\n",
                                                  " b: ", cd.gravity_wave_coefficient)

end # module