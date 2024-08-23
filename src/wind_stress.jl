module WindStressModel

export WindStress, WindStressBoundaryConditions, LogarithmicNeutralWind

using Roots

using Adapt: adapt

using Oceananigans.Architectures: on_architecture, CPU
using Oceananigans.BoundaryConditions: FluxBoundaryCondition
using Oceananigans.BuoyancyModels: g_Earth

using Walrus: get_value, normalise_surface_function
using Walrus.Interpolations: SimpleInterpolation, Limited

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

    relative_speed = wind_speed#√((uʷ - u)^2 + (vʷ - v)^2)

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

struct LogarithmicNeutralWind{FT, Z} # can't think of a good name for this
  monin_obukhov_stability_length :: FT
            charnock_coefficient :: FT
         air_kinematic_viscosity :: FT
        gravity_wave_coefficient :: FT
            gravity_acceleration :: FT

                roughness_length :: Z
end

"""
    LogarithmicNeutralWind(; monin_obukhov_stability_length = 0.4
                             charnock_coefficient = 0.014
                             air_kinematic_viscosity = 1.488e-5
                             gravity_wave_coefficient = 0.11
                             gravity = g_Earth,

                             precomputed_roughness_length = false,
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

                                  precomputed_roughness_length = false,
                                  precompute_wind_speeds = [0:25/100000:25;],
                                  arch = CPU()) where FT

    if precomputed_roughness_length
        tmp = LogarithmicNeutralWind(monin_obukhov_stability_length, charnock_coefficient,
                                     air_kinematic_viscosity, gravity_wave_coefficient, gravity_acceleration,
                                     nothing)

        κ = monin_obukhov_stability_length
        ν = air_kinematic_viscosity
        aᶜ = charnock_coefficient
        b = gravity_wave_coefficient
        g = gravity_acceleration

        params = (; κ, ν, aᶜ, b, g)

        n = 100000
        lengths = zeros(length(precompute_wind_speeds))

        for (n, wind_speed) in enumerate(precompute_wind_speeds)
            lengths[n] = find_velocity_roughness_length(tmp, wind_speed, 10, params)
        end

        roughness_length = SimpleInterpolation(precompute_wind_speeds, lengths; arch, mode = Limited())
    else
        roughness_length = nothing
    end

    return LogarithmicNeutralWind(monin_obukhov_stability_length, charnock_coefficient,
                                  air_kinematic_viscosity, gravity_wave_coefficient, gravity_acceleration,
                                  roughness_length)
end

adapt_structure(to, dc::LogarithmicNeutralWind) = 
    LogarithmicNeutralWind(dc.monin_obukhov_stability_length, dc.charnock_coefficient,
                           dc.air_kinematic_viscosity, dc.gravity_wave_coefficient, dc.gravity_acceleration,
                           adapt(to, dc.roughness_length))

@inline velocity_roughness_length_roots(z₀, params) = 
    log(params.reference_height/z₀)/(params.κ * params.wind_speed) * 
        (params.ν * params.b + params.aᶜ * (params.κ * params.wind_speed / log(params.reference_height/z₀))^3) - z₀

"""
    find_velocity_roughness_length(wind_speed, reference_height, params)

A function that finds the velocity roughness length for the `LogarithmicNeutralWind` drag
coefficient model.

This will sometimes fail as the function is not well behaved at either low reference heights
(it has been tuned for 10m wind), or high (⪆ 20 m/s).
"""
@inline function find_velocity_roughness_length(::LogarithmicNeutralWind{<:Any, Nothing}, wind_speed, reference_height, params)
    z₀ = reference_height
    
    upper_bounds_guess = ifelse(wind_speed < 0.05, 0.95 * reference_height, ifelse(wind_speed < 14, 0.5 * reference_height, 0.2 * reference_height))

    # do we want to tolarate this branch? Does it make it faster in the normal case?
    wind_speed < 0.01 || (z₀ = find_zero(velocity_roughness_length_roots, (0.00001, upper_bounds_guess), Bisection(), p = merge(params, (; wind_speed, reference_height))))

    return z₀
end

@inline find_velocity_roughness_length(dc::LogarithmicNeutralWind, wind_speed, reference_height, params) = dc.roughness_length(wind_speed)

@inline function (dc::LogarithmicNeutralWind)(wind_speed)
    κ = dc.monin_obukhov_stability_length
    ν = dc.air_kinematic_viscosity
    aᶜ = dc.charnock_coefficient
    b = dc.gravity_wave_coefficient
    g = dc.gravity_acceleration

    params = (; κ, ν, aᶜ, b, g)

    z₀ = find_velocity_roughness_length(dc, wind_speed, 10, params)

    Cᵈ = (params.κ / log(10 / z₀)) ^ 2

    Cᵈ = ifelse(isfinite(Cᵈ), Cᵈ, 0) # occurs when z₀ -> reference_height
    
    return Cᵈ
end

summary(::LogarithmicNeutralWind) = string("Log neutral wind drag coefficient model")
show(io::IO, cd::LogarithmicNeutralWind) = println(io, summary(cd), " with:\n",
                                                  " κ: ",  cd.monin_obukhov_stability_length, "\n",
                                                  " ν: ",  cd.air_kinematic_viscosity, " m²/s\n",
                                                  " aᶜ: ", cd.charnock_coefficient, "\n",
                                                  " g: " , cd.gravity_acceleration, "m/s²\n",
                                                  " b: ", cd.gravity_wave_coefficient)

end # module