using Oceananigans: fields
using Oceananigans.Architectures: architecture
using Oceananigans.BuoyancyFormulations: g_Earth
using Oceananigans.Fields: Field, Center, set!
using Oceananigans.Utils: launch!
using KernelAbstractions: @kernel, @index

struct SimilarityTheoryInterface{FT, VT, VP, SF, RL, DC, HC, MI} # can't think of a good name for this
             von_karman_constant :: FT
            gravity_acceleration :: FT
                reference_height :: FT
             virtual_temperature :: VT
   virtual_potential_temperature :: VP
           stability_formulation :: SF
                roughness_length :: RL

                drag_coefficient :: DC
       heat_exchange_coefficient :: HC

                  max_iterations :: MI
end

"""

"""
function SimilarityTheoryInterface(grid; 
                                   von_karman_constant::FT = 0.4,
                                   gravity_acceleration::FT = g_Earth,
                                   reference_height::FT = 10.0,
                                   virtual_temperature = VirtualTemperature(),
                                   virtual_potential_temperature = VirtualPotentialTemperature(),
                                   stability_parameterisation = DyerPaulsonStabilityFormulation(),
                                   roughness_length = SmoothAndCharnock(),
                                   max_iterations = 40) where FT

    drag_coefficient = Field{Center, Center, Nothing}(grid; indices = (:, :, 1))
    heat_exchange_coefficient = Field{Center, Center, Nothing}(grid; indices = (:, :, 1))

    set!(drag_coefficient, sqrt(1e-3))
    set!(heat_exchange_coefficient, sqrt(1e-3))

    return SimilarityTheoryInterface(von_karman_constant, gravity_acceleration, reference_height,
                                     virtual_temperature, virtual_potential_temperature,
                                     stability_parameterisation, roughness_length,
                                     drag_coefficient, heat_exchange_coefficient, 
                                     max_iterations)
end

adapt_structure(to, dc::SimilarityTheoryInterface) = 
    SimilarityTheoryInterface(dc.von_karman_constant, dc.gravity_acceleration, dc.reference_height,
                              adapt(to, dc.virtual_temperature), adapt(to, dc.virtual_potential_temperature),
                              adapt(to, dc.stability_formulation), adapt(to, dc.roughness_length),
                              adapt(to, dc.drag_coefficient), adapt(to, dc.heat_exchange_coefficient),
                              adapt(to, dc.max_iterations))

@inline function itterate_scaling_values(u′, T′, U, θ, T, w, zᵤ, zₜ, p)
    FT = typeof(T)

    Tᵥ = p.virtual_temperature(T + FT(273.15), w)

    θᵥ = (θ + FT(273.15)) * (1 + 4.7e-4 * w * (θ + FT(273.15)) / 273.15) + zₜ * 0.01

    κ = p.von_karman_constant
    g = p.gravity_acceleration

    Cₕ = -u′ * T′ / (U * (T - θ))

    L = -u′^3 * θᵥ / (g * κ * Cₕ * U * (Tᵥ - θᵥ))

    #L = ifelse(isinf(Cₕ), Inf, 0)

    zₒ, zₒₜ = p.roughness_length(abs(u′))

    ψₘ, _ = p.stability_formulation(zᵤ, L)
    _, ψₜ = p.stability_formulation(zₜ, L)

    ψₘ = ifelse(isnan(ψₘ), 0, ψₘ)
    ψₜ = ifelse(isnan(ψₜ), 0, ψₜ)

    u′₊ = κ * U / (log(zᵤ/zₒ) - ψₘ)
    T′₊ = κ * (θ - T) / (log(zₜ/zₒₜ) - ψₜ)

    # SMITH, 1988 says κ * (θᵥ - Tᵥ) / (log(zₜ/zₒₜ) - ψₜ) but that makes 
    # very strange values around θ -(⨥/−)→ T

    return (; u′ = u′₊, T′ = T′₊)
end

@kernel function _compute_coefficients!(interface::SimilarityTheoryInterface, grid, clock, model_fields, atmosphere)
    i, j = @index(Global, NTuple)

    FT = eltype(grid)

    U = relative_wind_speed(atmosphere, i, j, grid, clock, model_fields)
    θ = temperature(atmosphere, i, j, grid, clock, model_fields)
    w = air_water_mixing_ratio(atmosphere, i, j, grid, clock, model_fields)
    zᵤ = velocity_reference_height(atmosphere, i, j, grid, clock, model_fields)
    zₜ = temperature_reference_height(atmosphere, i, j, grid, clock, model_fields)
    T = @inbounds model_fields.T[i, j, grid.Nz]

    u′, T′ = sqrt(FT(1e-3)), sqrt(FT(1e-3))

    u′₋, T′₋ = FT(Inf), FT(Inf)

    iters = 0
    
    @inbounds while ((abs(u′ - u′₋) > 1e-8) | (abs(T′ - T′₋) > 1e-8)) & (iters <= interface.max_iterations)
        u′₋ = u′
        T′₋ = T′

        next_step = itterate_scaling_values(u′, T′, U, θ, T, w, zᵤ, zₜ, interface)

        u′ = next_step.u′
        T′ = next_step.T′

        iters += 1
    end

#    ((abs(u′ - u′₋) > 1e-8) | (abs(T′ - T′₋) > 1e-8)) && @warn "Did not converge with $u′, $T′, $U, $θ, $T, $w"

    Cd = @inbounds u′^2 / (U^2 + eps(0.0))
    Ch = @inbounds - T′ * u′ / (T - θ + eps(0.0)) / (U + eps(0.0))

    @inbounds interface.drag_coefficient[i, j, 1] = min(convert(FT, 1/10), ifelse(U == 0, zero(FT), Cd))
    @inbounds interface.heat_exchange_coefficient[i, j, 1] = min(convert(FT, 1/10), ifelse(isfinite(Ch), Ch, FT(1e-3)))
end

@inline function update_interface!(interface, model, atmosphere)
    clock = model.clock
    grid = model.grid
    model_fields = fields(model)
    arch = architecture(grid)

    launch!(arch, grid, :xy, _compute_coefficients!, interface, grid, clock, model_fields, atmosphere)

    return nothing
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
                       adapt(to, vpt.adiabatic_gradient))

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
@inline drag_coefficient(interface::SimilarityTheoryInterface, i, j, grid, clock, model_fields, atmosphere) =
    @inbounds interface.drag_coefficient[i, j, 1]

@inline heat_exchange_coefficient(interface::SimilarityTheoryInterface, i, j, grid, clock, model_fields, atmosphere) =
    @inbounds interface.heat_exchange_coefficient[i, j, 1]
