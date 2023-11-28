"""
TidalForcing

Provides quick setup of tidal forcing
"""
module TidalForcings

export Tide, TidalForcing

using Oceananigans.Coriolis: fᶠᶠᵃ, AbstractRotation
using Oceananigans.Forcings: Forcing
using Oceananigans.Units: hours

"""
    Tide(; x_amplitude, 
           y_amplitude, 
           period = 12.3782216453hours,
           nodal_time = 0., 
           x_lag = 0., 
           y_lag = 0.,
           coriolis = nothing)

Sets up a model of tidal forcing with default parameters of an M2 tide.

Keyword Arguments
=================
- `x_amplitude`: the tidal amplitude in the x direction
- `y_amplitude`: the tidal amplitude in the x direction
- `period`: the tidal period (defaults to that of an M2 tide)
- `nodal_time`: the time at which peak flow occurs
- `x_lag`: the phase lag for the tidal component in the x direction
- `y_lag`: the phase lag for the tidal component in the y direction
- `coriolis`: a model for the coriolis parameter 

Example
=======

```jldoctest
julia> using Walrus: Tide

julia> using Oceananigans

julia> tide = Tide(x_amplitude = 0.1, y_amplitude = 0.)
(::Tide{Float64, Nothing}) (generic function with 2 methods)
julia> forcing = (u = Forcing(tide, parameters = Val(:x), discrete_form = true),
                  v = Forcing(tide, parameters = Val(:y), discrete_form = true))
(u = DiscreteForcing{Val{:x}}
├── func: (::Tide{Float64, Nothing}) (generic function with 2 methods)
└── parameters: Val{:x}(), v = DiscreteForcing{Val{:y}}
├── func: (::Tide{Float64, Nothing}) (generic function with 2 methods)
└── parameters: Val{:y}())
```
"""
@kwdef struct Tide{FT, C} <: Function
    x_amplitude :: FT
    y_amplitude :: FT
         period :: FT = 12.3782216453hours
     nodal_time :: FT = 0.
          x_lag :: FT = 0.
          y_lag :: FT = 0.
 coriolis_model :: C  = nothing
end

adapt_structure(to, t::Tide) = Tide(t.x_amplitude,
                                    t.y_amplitude,
                                    t.period,
                                    t.nodal_time,
                                    t.x_lag,
                                    t.y_lag,
                                    adapt(to, t.coriolis))

function (tide::Tide)(i, j, k, grid, clock, model_fields, ::Val{:x})
    f = fᶠᶠᵃ(i, j, k, grid, tide.coriolis)

    ω = 2π / tide.period
    t = clock.time

    return - (tide.x_amplitude * ω * sin(ω * (t - tide.nodal_time) - tide.x_lag)
              + tide.y_amplitude * f * sin(ω * (t - tide.nodal_time) - tide.y_lag))
end

function (tide::Tide)(i, j, k, grid, clock, model_fields, ::Val{:y})
    f = fᶠᶠᵃ(i, j, k, grid, tide.coriolis)

    ω = 2π / tide.period
    t = clock.time

    return - (tide.y_amplitude * ω * sin(ω * (t - tide.nodal_time) - tide.y_lag)
              - tide.x_amplitude * f * sin(ω * (t - tide.nodal_time) - tide.x_lag))
end

"""
    TidalForcing(; x_amplitude,
                   y_amplitude,
                   period = 12.3782216453hours,
                   nodal_time = 0.,
                   x_lag = 0.,
                   y_lag = 0.,
                   coriolis = nothing)

A convenience constructor for `Tide` which returns the forcings pre wrapped.


Keyword Arguments
=================
- `x_amplitude`: the tidal amplitude in the x direction
- `y_amplitude`: the tidal amplitude in the x direction
- `period`: the tidal period (defaults to that of an M2 tide)
- `nodal_time`: the time at which peak flow occurs
- `x_lag`: the phase lag for the tidal component in the x direction
- `y_lag`: the phase lag for the tidal component in the y direction
- `coriolis`: a model for the coriolis parameter 

Example
=======

```jldoctest
julia> using Walrus: TidalForcing

julia> tidal_forcing = TidalForcing(x_amplitude = 0.1, y_amplitude = 0.)
(u = DiscreteForcing{Val{:x}}
├── func: (::Tide{Float64, Nothing}) (generic function with 2 methods)
└── parameters: Val{:x}(), v = DiscreteForcing{Val{:y}}
├── func: (::Tide{Float64, Nothing}) (generic function with 2 methods)
└── parameters: Val{:y}())
```
"""
function TidalForcing(; x_amplitude,
                        y_amplitude,
                        period = 12.3782216453hours,
                        nodal_time = 0.,
                        x_lag = 0.,
                        y_lag = 0.,
                        coriolis = nothing)

    tide = Tide(x_amplitude, y_amplitude,
                period, nodal_time,
                x_lag, y_lag, coriolis)

    return (u = Forcing(tide, parameters = Val(:x), discrete_form = true),
            v = Forcing(tide, parameters = Val(:y), discrete_form = true))
end
end # module