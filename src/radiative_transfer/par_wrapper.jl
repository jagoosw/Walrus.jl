# extends the OceanBioME light attenuation model to heat the water
# assumes that all absorbed light (including that absorbed by phytoplankton etc.) warms the water

using Oceananigans.Operators: ∂zᶜᶜᶜ

struct PARModelHeating{L, FT}
 light_attenuation_model :: L
     water_heat_capacity :: FT
           water_density :: FT
end

function PARModelHeating(; light_attenuation_model::L,
                           water_heat_capacity::FT = 3991.0, 
                           water_density::FT = 1026.0) where {L, FT}
    return PARModelHeating{L, FT}(light_attenuation_model, water_heat_capacity, water_density)
end

adapt_structure(to, bh::PARModelHeating) = 
    PARModelHeating(adapt(to, bh.light_attenuation_model),
                    bh.water_heat_capacity,
                    bh.water_density)
                    
@inline function (heating::PARModelHeating)(i, j, k, grid, clock, model_fields)
    ρ = heating.water_density
    cₚ = heating.water_heat_capacity

    light_attenuation = ∂zᶜᶜᶜ(i, j, k, grid, heating.field)

    return light_attenuation / (ρ * cₚ)
end
