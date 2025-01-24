module WindStressModel

export WindStress, WindStressBoundaryConditions, LogarithmicNeutralWind

using Roots

using Adapt: adapt

using Oceananigans.Architectures: on_architecture, CPU
using Oceananigans.BoundaryConditions: FluxBoundaryCondition
using Oceananigans.BuoyancyModels: g_Earth

using Walrus: get_value, normalise_surface_function
using Walrus.Interpolations: SimpleInterpolation

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
                      drag_coefficient = LogarithmicNeutralWind(), 
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

    stress_velocity = ρₐ / ρₒ * wind_stress.drag_coefficient(relative_speed) * relative_speed

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

struct LogarithmicNeutralWind{FT, CD} # can't think of a good name for this
  monin_obukhov_stability_length :: FT
            charnock_coefficient :: FT
         air_kinematic_viscosity :: FT
        gravity_wave_coefficient :: FT
            gravity_acceleration :: FT

                drag_coefficient :: CD
end

"""
    LogarithmicNeutralWind(; monin_obukhov_stability_length = 0.4
                             charnock_coefficient = 0.014
                             air_kinematic_viscosity = 1.488e-5
                             gravity_wave_coefficient = 0.11
                             gravity = g_Earth,

                             precompute_drag_coefficients = false,
                             precompute_wind_speeds = [0:25/100000:25;],
                             arch = CPU())

Returns a `LogarithmicNeutralWind` parameterisation for the surface drag coefficient

``C_d`` is parameterised as,
```math
C_d = \\left(\\frac{\\kappa}{\\log{\\frac{10}{z_0}}}\\right)^2,
```
where ``\\kappa`` is the Monin‐Obukhov stability length and ``z_0`` is the velocity 
roughness length. This is the roughness length scale which logarithmically brings 
the relative velocity to zero at the surface, i.e.
```math
U=\\frac{u\\star}{\\kappa}\\log\\frac{z}{z_0},
```
where ``u\\star`` is the friction velocity. Additionally ``z_0`` is given as,
```math
z_0=b\\frac{\\nu}{u\\star} + \\frac{a_c}{g}u\\star^2,
```
where ``\\nu`` is the kinematic viscosity of air and g is the acceleration of gravity.

This model itterativly solves these equations to find ``z_0``. Alternativly, if the flag 
`precomputed_roughness_length` is set to they are pre computed at `precompute_wind_speeds` 
between which ``z_0`` is then interpolated during run time. Precomputed velocities are 
converted to appropriate types for `arch` (i.e. `CPU()` or `GPU()`)

This parameterisaion is described in [smith1988](@citet)
"""
function LogarithmicNeutralWind(; monin_obukhov_stability_length::FT = 0.4,
                                  charnock_coefficient::FT = 0.014,
                                  air_kinematic_viscosity::FT = 1.488e-5,
                                  gravity_wave_coefficient::FT = 0.11,
                                  gravity_acceleration::FT = g_Earth,

                                  precompute_drag_coefficients = false,
                                  precompute_wind_speeds = [0:1e-4:100;],
                                  arch = CPU()) where FT

    if precompute_drag_coefficients
        tmp_Cd = LogarithmicNeutralWind(monin_obukhov_stability_length, charnock_coefficient,
                                        air_kinematic_viscosity, gravity_wave_coefficient, gravity_acceleration,
                                        nothing)


        Cd = tmp_Cd.(precompute_wind_speeds)

        drag_coefficient = SimpleInterpolation(precompute_wind_speeds, Cd; arch)
    else
        drag_coefficient = nothing
    end

    return LogarithmicNeutralWind(monin_obukhov_stability_length, charnock_coefficient,
                                  air_kinematic_viscosity, gravity_wave_coefficient, gravity_acceleration,
                                  drag_coefficient)
end

adapt_structure(to, dc::LogarithmicNeutralWind) = 
    LogarithmicNeutralWind(dc.monin_obukhov_stability_length, dc.charnock_coefficient,
                           dc.air_kinematic_viscosity, dc.gravity_wave_coefficient, dc.gravity_acceleration,
                           adapt(to, dc.drag_coefficient))

@inline function drag_coefficient_excess(Cd, p)
    κ, ν, α, b, g, u = p

    u = max(0.01, u)

    u′ = √(Cd * u^2)

    z₀ = b * ν / u′ + α * u′^2 / g

    Cd_new = (κ / log(10 / z₀))^2

    return Cd_new - Cd
end

@inline (dc::LogarithmicNeutralWind)(wind_speed) = dc.drag_coefficient(wind_speed)

@inline function (dc::LogarithmicNeutralWind{<:Any, Nothing})(wind_speed)
    κ = dc.monin_obukhov_stability_length
    ν = dc.air_kinematic_viscosity
    α = dc.charnock_coefficient
    b = dc.gravity_wave_coefficient
    g = dc.gravity_acceleration

    p = (; κ, ν, α, b, g, wind_speed)
    
    return find_zero(drag_coefficient_excess, (0.0005, 0.02), Bisection(), p)
end

summary(::LogarithmicNeutralWind) = string("Log neutral wind drag coefficient model")
show(io::IO, cd::LogarithmicNeutralWind) = println(io, summary(cd), " with:\n",
                                                  " κ: ",  cd.monin_obukhov_stability_length, "\n",
                                                  " ν: ",  cd.air_kinematic_viscosity, " m²/s\n",
                                                  " aᶜ: ", cd.charnock_coefficient, "\n",
                                                  " g: " , cd.gravity_acceleration, "m/s²\n",
                                                  " b: ", cd.gravity_wave_coefficient)

end # module
