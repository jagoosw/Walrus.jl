"""
WindDrivenStokesParameterisation

A stokes drift parameterisation based on wind speed where the peak frequency is determined
from [PiersonMoskowitz](@citet) where:
```math
\\omega_p=\\frac{2\\pi g}{7.69 u_{10}}.
```

The significant wave height is then determined from "Toba's 3/2 power law" [Toba](@citep):
```math
H = 0.057\\frac{4u^\\star^2}{g}\\left(\\frac{u^\\star\\omega_p}{g}\\right)^{-3/2}.
```

And finally the depth attenuation is determined from as:
```math
\\omega = \\sqrt{gk\\tanh hk},
```
which is approximated to:
```math
\\omega = \\sqrt{gk},
```
for `h * k_0 > π` where `\\omega = \\sqrt{gk_0}`.

Given that the peak frequency is for fully developed sea states this parameterisation is only valid
for cases where the wind speed evolves slowly compared to the wave age. Additionally, this 
parameterisation is only for monochromatic wave fields. Since stokes drift in Oceananigans doesn't
have access to the velocity fields this model also has to neglect the developement of the sea state
with respect to the wind (i.e. can't work out the relative wind speed and instead has to just use
the wind speed).
"""    

module WindDrivenStokesParameterisation

export WindDrivenStokesDrift, WindDrivenStokesDriftSetup

using Roots

using Oceananigans.Architectures: on_architecture, CPU, architecture
using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.StokesDrifts: UniformStokesDrift

using Walrus.Interpolations: SimpleInterpolation

import Adapt: adapt_structure
import Base: summary, show

struct WindDrivenStokesDrift{DI, WI, DE, WN, G, TW}
                   direction :: DI
                        wind :: WI
                       depth :: DE
                 wavenumbers :: WN
  gravitational_acceleration :: G
   time_interpolation_window :: TW
end

function WindDrivenStokesDrift(; wind, depth,
                                 direction = 0,
                                 gravitational_acceleration = g_Earth,
                                 time_interpolation_window = 240,
                                 precomputed_wavenumbers = false,
                                 precomputed_peak_frequencies = [0.3:0.001:1000;],
                                 grid = nothing,
                                 arch = isnothing(grid) ? CPU() : architecture(grid))
    if precomputed_wavenumbers
        precomputed_k = similar(precomputed_peak_frequencies)
        g = gravitational_acceleration
        h = depth

        for (n, ωₚ) in enumerate(precomputed_peak_frequencies)
            precomputed_k[n] = find_zero(wavenumber, (0, 8 * ωₚ ^ 2 / g), Bisection(), p = (; g, ωₚ, h))
        end

        precomputed_wavenumbers = SimpleInterpolation(precomputed_peak_frequencies, precomputed_k; arch)
    end

    return WindDrivenStokesDrift(direction, wind, depth, precomputed_wavenumbers, gravitational_acceleration, time_interpolation_window)
end

adapt_structure(to, sd::WindDrivenStokesDrift) = 
    WindDrivenStokesDrift(adapt(to, sd.direction),
                          adapt(to, sd.wind),
                          adapt(to, sd.depth),
                          adapt(to, sd.wavenumbers),
                          adapt(to, sd.gravitational_acceleration),
                          adapt(to, sd.time_interpolation_window))

@inline shallow_water_attenuation(z, params) = ifelse(params.wavenumber * params.depth > 4, 
                                                      exp(params.wavenumber * z),
                                                      cosh(params.wavenumber * (z + params.depth)) / sinh(params.wavenumber * params.depth))

@inline ∂z_shallow_water_attenuation(z, params) = ifelse(params.wavenumber * params.depth > 4, 
                                                         params.wavenumber * exp(params.wavenumber * z),
                                                         params.wavenumber * sinh(params.wavenumber * (z + params.depth)) / sinh(params.wavenumber * params.depth))

@inline wavenumber(k, params) = params.ωₚ - √(params.g * k * tanh(params.h * k))

@inline find_wavenumber(sd::WindDrivenStokesDrift, ωₚ, k₀) = sd.wavenumbers(ωₚ)

@inline function find_wavenumber(sd::WindDrivenStokesDrift{<:Any, <:Any, <:Any, Nothing}, ωₚ, k₀)
    g = sd.gravitational_acceleration
    h = sd.depth

    return find_zero(wavenumber, (0, 8 * ωₚ ^ 2 / g), Bisection(), p = (; g, ωₚ, h))
end

@inline function surface_drift_velocity(z, t, params)
    wind = params.wind

    uʷ = abs(wind.reference_wind_speed(0, 0, t))

    dc = wind.drag_coefficient

    ρₐ = wind.air_density
    ρₒ = wind.water_density

    uₛ = sqrt(ρₐ / ρₒ * dc(uʷ) * uʷ ^ 2)

    g = params.gravitational_acceleration
    h = params.depth

    ωₚ = 2π * g / (7.69 * uʷ)

    H = 4 * uₛ ^ 2 / g * 0.057 * (uₛ * ωₚ / g) ^ (-3/2)

    k₀ = ωₚ ^ 2 / g

    (2π / k₀ < 2 * h) ?
        k = k₀ : 
        k = find_wavenumber(params, ωₚ, k₀)

    uₛ = H * ωₚ

    return uₛ, k
end

@inline function drift_velocity(z, t, params)
    uₛ, wavenumber = surface_drift_velocity(z, t, params)
    
    return uₛ * shallow_water_attenuation(z, (; depth = params.depth, wavenumber))
end

@inline function ∂z_drift_velocity(z, t, params)
    uₛ, wavenumber = surface_drift_velocity(z, t, params)
    
    return uₛ * ∂z_shallow_water_attenuation(z, (; depth = params.depth, wavenumber))
end

@inline uˢ(z, t, params) = - drift_velocity(z, t, params) * sind(params.direction)
@inline ∂t_uˢ(z, t, params) = (uˢ(z, t + params.time_interpolation_window/2, params) - uˢ(z, t - params.time_interpolation_window/2, params))/params.time_interpolation_window
@inline ∂z_uˢ(z, t, params) = - ∂z_drift_velocity(z, t, params) * sind(params.direction)

@inline vˢ(z, t, params) = - drift_velocity(z, t, params) * cosd(params.direction)
@inline ∂t_vˢ(z, t, params) = (vˢ(z, t + params.time_interpolation_window/2, params) - vˢ(z, t - params.time_interpolation_window/2, params))/params.time_interpolation_window
@inline ∂z_vˢ(z, t, params) = - ∂z_drift_velocity(z, t, params) * cosd(params.direction)

function WindDrivenStokesDriftSetup(; wind, depth,
                                        direction = 0,
                                        gravitational_acceleration = g_Earth,
                                        time_interpolation_window = 240,
                                        precomputed_wavenumbers = false,
                                        precomputed_peak_frequencies = [0.3:0.001:1000;],
                                        grid = nothing,
                                        arch = isnothing(grid) ? CPU() : architecture(grid))

    parameterisation = WindDrivenStokesDrift(; wind, depth, direction, gravitational_acceleration, time_interpolation_window, precomputed_wavenumbers, precomputed_peak_frequencies, grid, arch)

    return  UniformStokesDrift(; ∂z_uˢ, ∂z_vˢ, ∂t_uˢ, ∂t_vˢ, parameters = parameterisation)
end

end # module