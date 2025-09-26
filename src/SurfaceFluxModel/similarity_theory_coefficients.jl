using Oceananigans.BuoyancyFormulations: g_Earth

struct SimilarityTheoryInterface{FT, VT, VP, SF, RL} # can't think of a good name for this
             von_karman_constant :: FT
            gravity_acceleration :: FT
                reference_height :: FT
             virtual_temperature :: VT
   virtual_potential_temperature :: VP
           stability_formulation :: SF
                roughness_length :: RL
end

"""

"""
function SimilarityTheoryInterface(; von_karman_constant::FT = 0.4,
                                  gravity_acceleration::FT = g_Earth,
                                  reference_height::FT = 10.0,
                                  virtual_temperature = VirtualTemperature(),
                                  virtual_potential_temperature = VirtualPotentialTemperature(),
                                  stability_parameterisation = DyerPaulsonStabilityFormulation(),
                                  roughness_length = SmoothAndCharnock()) where FT

    return SimilarityTheoryInterface(von_karman_constant, gravity_acceleration, reference_height,
                                     virtual_temperature, virtual_potential_temperature,
                                     stability_parameterisation, roughness_length)
end

adapt_structure(to, dc::SimilarityTheoryInterface) = 
    SimilarityTheoryInterface(dc.von_karman_constant, dc.gravity_acceleration, dc.reference_height,
                              adapt(to, dc.virtual_temperature), adapt(to, dc.virtual_potential_temperature),
                              adapt(to, dc.stability_parameterisation), adapt(to, dc.roughness_length))

@inline function itterate_scaling_values(interface, previous_values, p)
    FT = typeof(interface.T)

    u′, T′ = previous_values

    T  = interface.T
    w  = interface.w
    Tᵥ = p.virtual_temperature(T + FT(273.15), w)

    θ  = interface.θ
    θᵥ = p.virtual_potential_temperature(θ + FT(273.15), w)

    U = interface.U
    κ = p.von_karman_constant
    g = p.gravity_acceleration
    zᵤ = interface.zᵤ
    zₜ = interface.zₜ

    Cₕ = -u′ * T′ / (U * (T - θ))

    Cₕ = ifelse(isinf(Cₕ), FT(1e-3), Cₕ)

    L = -u′^3 * Tᵥ / (g * κ * Cₕ * U * (Tᵥ - θᵥ))

    L = ifelse(isinf(Cₕ) | isnan(L), zero(T), L)

    zₒ, zₒₜ = p.roughness_length(u′)

    ψₘ, _ = p.stability_formulation(zᵤ, L)
    _, ψₜ = p.stability_formulation(zₜ, L)
    ψₘₒ, ψₜₒ = p.stability_formulation(zₒ, L)
    ψₘₒ, ψₜₒ = p.stability_formulation(zₒₜ, L)

    u′₊ = κ * U / (log(zᵤ/zₒ) - ψₘ + ψₘₒ)
    T′₊ = κ * (θᵥ - Tᵥ) / (log(zₜ/zₒₜ) - ψₜ + ψₜₒ)

    u′₊ = ifelse(isinf(zₒ), 0, u′₊)
    T′₊ = ifelse(isinf(zₒₜ)|isinf(zₒ), 0, T′₊)

    return max(0, u′₊), T′₊
end

@inline function (cd::SimilarityTheoryInterface)(U, θ, T, w, zᵤ, zₜ)
    interface = (; U, θ, T, w, zᵤ, zₜ)
    
    u′, T′ = sqrt.(1e-3), sqrt.(1e-3)

    u′₋, T′₋ = Inf, Inf

    iters = 0

    while ((abs(u′ - u′₋) > 1e-8) | (abs(T′ - T′₋) > 1e-8)) & (iters <= 20)
        u′₋ = u′
        T′₋ = T′

        u′, T′ = itterate_scaling_values(interface, (; u′, T′), cd)

        iters += 1
    end

    ((abs(u′ - u′₋) > 1e-8) | (abs(T′ - T′₋) > 1e-8)) && @warn "$cd coefficients did not converge"

    Cd = u′^2 / (U^2 + eps(0.0))
    Ch = - T′ * u′ / (T - θ + eps(0.0)) / (U + eps(0.0))
    
    Cd = ifelse(U == 0, 0, Cd)
    Ch = ifelse(T == θ, 0, Ch)

    return Cd, Ch
end

summary(::SimilarityTheoryInterface) = string("Similarity theory drag and heat transfer coefficients")
#=show(io::IO, cd::SimilarityTheoryInterface) = println(io, summary(cd), " with:\n",
                                                  " κ: ",  cd.von_karman_constant, "\n",
                                                  " g: " , cd.gravity_acceleration, "m/s²\n",
                                                  " ψ: ", cd.stability_formulation, "\n",
                                                  " zₒ: ", cd.rough)
=#
@kwdef struct VirtualTemperature{AM, MR}
    air_mixing_ratio :: AM = 0.0102
          mass_ratio :: MR = 0.622
end

Adapt.adapt_structure(to, vt::VirtualTemperature) =
    VirtualTemperature(adapt(to, vt.air_mixing_ratio),
                       adapt(to, vt.mass_ratio))

@inline function (vt::VirtualTemperature)(T, w)
    ϵ = vt.mass_ratio

    return T * (w + ϵ) / (ϵ * (1 + w))
end

@inline (vt::VirtualTemperature)(T, ::Nothing) = vt(T, vt.air_mixing_ratio) # todo: should be get_value but cba getting all the info in here

@kwdef struct VirtualPotentialTemperature{FT}
                height :: FT = 10.0
    adiabatic_gradient :: FT = 0.01
end

Adapt.adapt_structure(to, vpt::VirtualPotentialTemperature) =
    VirtualTemperature(adapt(to, vpt.height),
                       adapt(tp, vpt.adiabatic_gradient))

@inline function (vpt::VirtualPotentialTemperature)(θ, w)
    z = vpt.height
    Γ = vpt.adiabatic_gradient

    return θ * (1 + 4.7e-4 * w * θ / 273.15) + z * Γ # query shouldn't this be -zΓ since its brought down from the reference height
end

@kwdef struct DyerPaulsonStabilityFormulation{FT}
      stable_coefficient :: FT = 5.0
    unstable_coefficient :: FT = 16.0
end

@inline function (sf::DyerPaulsonStabilityFormulation)(z, L)
    α = sf.stable_coefficient
    ψ⁺ = -α * z / L

    β = sf.unstable_coefficient
    x = (1-β*min(0, z/L))^(1/4) 

    ψ⁻ₘ = 2 * log((1+x)/2) + log((1+x^2)/2) - 2 * atan(x) + π/2 # Paulson 1970 
    ψ⁻ₜ = 2 * log((1+x^2)/2)

    return ifelse(L>0, (ψ⁺, ψ⁺), (ψ⁻ₘ, ψ⁻ₜ))
end

@kwdef struct SmoothAndCharnock{FT}
        charnock_coefficient :: FT = 0.014
     air_kinematic_viscosity :: FT = 14.88e-6
    gravity_wave_coefficient :: FT = 0.11
        gravity_acceleration :: FT = g_Earth
end

@inline function (z::SmoothAndCharnock)(u′)
    ν = z.air_kinematic_viscosity
    α = z.gravity_wave_coefficient
    a = z.charnock_coefficient
    g = z.gravity_acceleration

    zₒ = α * ν / u′ + a * u′^2 / g
    
    Re = u′ * zₒ / ν
    zₒₜ = min(1.15e-4, 5.5e-5 * Re ^ -0.6) # This is an empirical parameterisation by Fairall et al. 2003

    return zₒ, zₒₜ
end

#####
##### external interface
#####
@inline function drag_coefficient(interface::SimilarityTheoryInterface, i, j, grid, clock, model_fields, atmosphere)
    Cd, _ = coefficients(interface::SimilarityTheoryInterface, i, j, grid, clock, model_fields, atmosphere)

    return Cd
end

@inline function heat_exchange_coefficient(interface::SimilarityTheoryInterface, i, j, grid, clock, model_fields, atmosphere)
    _, Ch = coefficients(interface::SimilarityTheoryInterface, i, j, grid, clock, model_fields, atmosphere)

    return Ch
end

@inline function coefficients(interface::SimilarityTheoryInterface, i, j, grid, clock, model_fields, atmosphere)
    U = relative_wind_speed(atmosphere, i, j, grid, clock, model_fields)
    θ = temperature(atmosphere, i, j, grid, clock, model_fields)
    w = air_water_mixing_ratio(atmosphere, i, j, grid, clock, model_fields)
    zᵤ = velocity_reference_height(atmosphere, i, j, grid, clock, model_fields)
    zₜ = temperature_reference_height(atmosphere, i, j, grid, clock, model_fields)
    T = @inbounds model_fields.T[i, j, grid.Nz]
    
    return interface(U, θ, T, w, zᵤ, zₜ)
end