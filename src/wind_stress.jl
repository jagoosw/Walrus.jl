module WindStressModel

export WindStress, WindStressBoundaryConditions

using Roots

using Adapt: adapt

using Oceananigans.Architectures: on_architecture, CPU
using Oceananigans.BoundaryConditions: FluxBoundaryCondition
using Oceananigans.BuoyancyFormulations: g_Earth

using Walrus: get_value, normalise_surface_function
using Walrus.Interpolations: SimpleInterpolation
using Walrus.InterfaceCoefficients: SimilarityTheoryInterface

import Adapt: adapt_structure
import Base: summary, show

struct WindStress{WS, WD, DC, FT} <: Function
      reference_wind_speed :: WS
  reference_wind_direction :: WD
          drag_coefficient :: DC
               air_density :: FT
             water_density :: FT
end

adapt_structure(to, ws::WindStress) = WindStress(adapt(to, ws.reference_wind_speed),
                                                 adapt(to, ws.reference_wind_direction),
                                                 adapt(to, ws.drag_coefficient),
                                                 ws.air_density,
                                                 ws.water_density)

"""
    WindStress(; reference_wind_speed, 
                 reference_wind_direction,
                 drag_coefficient = LogarithmicNeutralWind(), 
                 air_density = 1.225, 
                 water_density = 1026.)

Returns a wind stress model where the stress is given by,
```math
\\frac{\\tau}{\\rho_o} = \\rho_aC_dSU_{x/y},
```
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
(::WindStress{Float64, Float64, LogarithmicNeutralWind{Float64, Nothing}, Float64}) (generic function with 2 methods)

julia> boundary_conditions = (u = FieldBoundaryConditions(top = FluxBoundaryCondition(wind_stress, parameters = Val(:x))),
                              v = FieldBoundaryConditions(top = FluxBoundaryCondition(wind_stress, parameters = Val(:y))))
(u = Oceananigans.FieldBoundaryConditions, with boundary conditions
├── west: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── east: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── south: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── north: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── bottom: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── top: FluxBoundaryCondition: ContinuousBoundaryFunction (::WindStress{Float64, Float64, LogarithmicNeutralWind{Float64, Nothing}, Float64}) at (Nothing, Nothing, Nothing)
└── immersed: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing), v = Oceananigans.FieldBoundaryConditions, with boundary conditions
├── west: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── east: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── south: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── north: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── bottom: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── top: FluxBoundaryCondition: ContinuousBoundaryFunction (::WindStress{Float64, Float64, LogarithmicNeutralWind{Float64, Nothing}, Float64}) at (Nothing, Nothing, Nothing)
└── immersed: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing))

```
"""
function WindStress(; reference_wind_speed, 
                      reference_wind_direction,
                      interface = SimilarityTheoryInterface(), 
                      air_density = 1.225, 
                      water_density = 1026.)
        
    reference_wind_speed = normalise_surface_function(reference_wind_speed)
    reference_wind_direction = normalise_surface_function(reference_wind_direction)

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

julia> using Oceananigans

julia> wind_stress_boundary_conditions = WindStressBoundaryConditions(; reference_wind_speed = 0.1, reference_wind_direction = 90.)
(u = FluxBoundaryCondition: DiscreteBoundaryFunction (::WindStress{Float64, Float64, LogarithmicNeutralWind{Float64, Nothing}, Float64}) with parameters Val{:x}, v = FluxBoundaryCondition: DiscreteBoundaryFunction (::WindStress{Float64, Float64, LogarithmicNeutralWind{Float64, Nothing}, Float64}) with parameters Val{:y})

julia> boundary_conditions = (u = FieldBoundaryConditions(top = wind_stress_boundary_conditions.u),
                              v = FieldBoundaryConditions(top = wind_stress_boundary_conditions.v))
(u = Oceananigans.FieldBoundaryConditions, with boundary conditions
├── west: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── east: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── south: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── north: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── bottom: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── top: FluxBoundaryCondition: DiscreteBoundaryFunction (::WindStress{Float64, Float64, LogarithmicNeutralWind{Float64, Nothing}, Float64}) with parameters Val{:x}
└── immersed: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing), v = Oceananigans.FieldBoundaryConditions, with boundary conditions
├── west: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── east: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── south: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── north: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── bottom: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── top: FluxBoundaryCondition: DiscreteBoundaryFunction (::WindStress{Float64, Float64, LogarithmicNeutralWind{Float64, Nothing}, Float64}) with parameters Val{:y}
└── immersed: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing))

```
"""
function WindStressBoundaryConditions(; reference_wind_speed, 
                                        reference_wind_direction,
                                        drag_coefficient = LogarithmicNeutralWind(), 
                                        air_density = 1.225, 
                                        water_density = 1026.)

    wind_stress = WindStress(; reference_wind_speed, 
                               reference_wind_direction,
                               drag_coefficient, air_density, water_density)

    u = FluxBoundaryCondition(wind_stress, parameters = Val(:x), discrete_form=true)

    v = FluxBoundaryCondition(wind_stress, parameters = Val(:y), discrete_form=true)

    return (; u, v)
end

@inline function (wind_stress::WindStress)(i, j, grid, clock, model_fields, ::Val{:x})
    ρₐ = wind_stress.air_density
    ρₒ = wind_stress.water_density

    t = clock.time

    u = @inbounds model_fields.u[i, j, grid.Nz]
    v = @inbounds model_fields.v[i, j, grid.Nz]

    wind_speed = get_value(wind_stress.reference_wind_speed, i, j, grid, clock)
    wind_direction = get_value(wind_stress.reference_wind_direction, i, j, grid, clock)

    uʷ = - wind_speed * sind(wind_direction)
    vʷ = - wind_speed * cosd(wind_direction)

    relative_speed = √((uʷ - u)^2 + (vʷ - v)^2)

    Cd, _ = wind_stress.drag_coefficient(relative_speed, )

    stress_velocity = ρₐ / ρₒ * Cd * relative_speed

    return - stress_velocity * (uʷ - u)
end

@inline function (wind_stress::WindStress)(i, j, grid, clock, model_fields, ::Val{:y})
    ρₐ = wind_stress.air_density
    ρₒ = wind_stress.water_density

    t = clock.time

    u = @inbounds model_fields.u[i, j, grid.Nz]
    v = @inbounds model_fields.v[i, j, grid.Nz]

    wind_speed = get_value(wind_stress.reference_wind_speed, i, j, grid, clock)
    wind_direction = get_value(wind_stress.reference_wind_direction, i, j, grid, clock)

    uʷ = - wind_speed * sind(wind_direction)
    vʷ = - wind_speed * cosd(wind_direction)

    relative_speed = √((uʷ - u)^2 + (vʷ - v)^2)

    stress_velocity = ρₐ / ρₒ * wind_stress.drag_coefficient(relative_speed) * relative_speed

    return - stress_velocity * (vʷ - v)
end

summary(::WindStress) = string("Wind stress model")
show(io::IO, wind::WindStress) = println(io, summary(wind), " with:\n",
                                     " Wind speed: ", summary(wind.reference_wind_speed), "\n",
                                     " Wind direction: ", summary(wind.reference_wind_direction), "\n",
                                     " Drag coefficient: ", summary(wind.drag_coefficient), "\n",
                                     " Air density: ", wind.air_density, " kg/m³\n",
                                     " Water density: ", wind.water_density, " kg/m³")

end # module