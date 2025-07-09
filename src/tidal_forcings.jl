module TidalForcing

export Tide, Tides, M2Tide

using Adapt
using Oceananigans.Forcings: Forcing
using Oceananigans.Coriolis: fᶠᶠᵃ, AbstractRotation
using Oceananigans.Units: hours

import Adapt: adapt_structure

struct Tide{FT, CO} <: Function
    x_amplitude :: FT 
    y_amplitude :: FT

       x_offset :: FT
       y_offset :: FT

         period :: FT

       coriolis :: CO
end

Adapt.adapt_structure(to, tide::Tide) =  
    Tide(adapt(to, tide.x_amplitude),
         adapt(to, tide.y_amplitude),
         adapt(to, tide.x_offset),
         adapt(to, tide.y_offset),
         adapt(to, tide.period),
         adapt(to, tide.coriolis))

function Tides(x_amplitude, y_amplitude, x_offset, y_offset, period, coriolis)
    tide = Tide(x_amplitude, y_amplitude, x_offset, y_offset, period, coriolis)

    x = Forcing(tide, discrete_form = true, parameters = Val(:x))
    y = Forcing(tide, discrete_form = true, parameters = Val(:y))

    return (; x, y)
end

function M2Tide(Ax, Ay; x_offset = 0.0, y_offset = 0.0, period = 12.4206012hours, coriolis)
    tide = Tide(Ax, Ay, x_offset, y_offset, period, coriolis)

    x = Forcing(tide, discrete_form = true, parameters = Val(:x))
    y = Forcing(tide, discrete_form = true, parameters = Val(:y))

    return (; x, y)
end

# single component

function (tide::Tide)(i, j, k, grid, clock, model_fields, ::Val{:x})
    f = fᶠᶠᵃ(i, j, k, grid, tide.coriolis)

    ω = 2π / tide.period
    t = clock.time

    return - (tide.x_amplitude * ω * sin(ω * t - tide.x_offset)
            + tide.y_amplitude * f * cos(ω * t - tide.y_offset))
end

function (tide::Tide)(i, j, k, grid, clock, model_fields, ::Val{:y})
    f = fᶠᶠᵃ(i, j, k, grid, tide.coriolis)

    ω = 2π / tide.period
    t = clock.time

    return - (tide.y_amplitude * ω * sin(ω * t - tide.y_offset)
            - tide.x_amplitude * f * cos(ω * t - tide.x_offset))
end

# 2 components
function (tide::Tide{<:NTuple{2}})(i, j, k, grid, clock, model_fields, ::Val{:x})
    f = fᶠᶠᵃ(i, j, k, grid, tide.coriolis)

    ω₁ = @inbounds 2π / tide.period[1]
    ω₂ = @inbounds 2π / tide.period[2]

    t = clock.time

    T₁ = @inbounds - (tide.x_amplitude[1] * ω₁ * sin(ω₁ * t - tide.x_offset[1]) + 
                      tide.y_amplitude[1] * f * cos(ω₁ * t - tide.y_offset[1]))
    T₂ = @inbounds - (tide.x_amplitude[2] * ω₂ * sin(ω₂ * t - tide.x_offset[2]) + 
                      tide.y_amplitude[2] * f * cos(ω₂ * t - tide.y_offset[2]))

    return T₁ + T₂
end

function (tide::Tide{<:NTuple{2}})(i, j, k, grid, clock, model_fields, ::Val{:y})
    f = fᶠᶠᵃ(i, j, k, grid, tide.coriolis)

    ω₁ = @inbounds 2π / tide.period[1]
    ω₂ = @inbounds 2π / tide.period[2]

    t = clock.time

    T₁ = @inbounds - (tide.y_amplitude[1] * ω₁ * sin(ω₁ * t - tide.y_offset[1]) - tide.x_amplitude[1] * f * cos(ω₁ * t - tide.x_offset[1]))
    T₂ = @inbounds - (tide.y_amplitude[2] * ω₂ * sin(ω₂ * t - tide.y_offset[2]) - tide.x_amplitude[2] * f * cos(ω₂ * t - tide.x_offset[2]))

    return T₁ + T₂
end

# n-component
@inline function x_tide_compoent(tide, t, f, n) 
    ω = @inbounds 2π / tide.period[n]

    return @inbounds - (tide.x_amplitude[n] * ω * sin(ω * t - tide.x_offset[n]) +
                        tide.y_amplitude[n] * ω * cos(ω * t - tide.y_offset[n]))
end

@inline function y_tide_compoent(tide, t, f, n) 
    ω = @inbounds 2π / tide.period[n]

    return @inbounds - (tide.y_amplitude[n] * ω * sin(ω * t - tide.y_offset[n]) -
                        tide.x_amplitude[n] * ω * cos(ω * t - tide.x_offset[n]))
end

function (tide::Tide{<:NTuple{N}})(i, j, k, grid, clock, model_fields, ::Val{:x}) where N
    f = fᶠᶠᵃ(i, j, k, grid, tide.coriolis)

    t = clock.time

    return sum(ntuple(n -> x_tide_compoent(tide, t, f, n), Val(N)))
end

function (tide::Tide{<:NTuple{N}})(i, j, k, grid, clock, model_fields, ::Val{:y}) where N
    f = fᶠᶠᵃ(i, j, k, grid, tide.coriolis)

    t = clock.time

    return sum(ntuple(n -> x_tide_compoent(tide, t, f, n), Val(N)))
end

end # module
