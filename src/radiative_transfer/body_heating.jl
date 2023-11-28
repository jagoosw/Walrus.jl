"""
    BodyHeating

A model for single band light attenuation which heats the water in the form:
``\\frac{\\partial I}{\\partial z} = -\\alpha I``
where ``I`` is the radiation intensity and ``\\alpha`` is the attenuation coefficient. This heats the water
like ``\\frac{\\partial T}{\\partial t} = \\frac{\\alpha I}{c^p\\rho}`` where ``c^p`` is the specific heat capacity
and ``\\rho`` is the water density. 
"""
struct BodyHeating{FT, S, F} <: AbstractBiogeochemistry
 water_attenuation_coefficient :: FT
           water_heat_capacity :: FT
                 water_density :: FT

                  surface_flux :: S
                         field :: F 
   
   BodyHeating(water_attenuation_coefficient::FT,
               water_heat_capacity::FT,
               water_density::FT,
               surface_flux::S,
               field::F) where {FT, S, F} = new{FT, S, F}(water_attenuation_coefficient,
                                                          water_heat_capacity,
                                                          water_density,
                                                          surface_flux,
                                                          field)
end

"""
    BodyHeating(; surface_flux,
                  grid = nothing,
                  field = nothing,
                  water_attenuation_coefficient = 1.8,
                  water_heat_capacity = 3991.0, # J K⁻¹ kg⁻¹
                  water_density = 1026.0) # kg m⁻³

Creates a model in which a `surface_flux` (W / m²) is attenuated by and heats the water.
This interacts with Oceananigans through the biogeochemistry but can also be used as
a modifier in `OceanBioME`'s `Biogeochemistry`.

Keyword Arguments
=================

- `surface_flux` (required): a function returning the surface radiaiton flux in the form `surface_flux(x, y, t)`
- `grid`: the model grid required to setup the radiaiton field (reccommended)
- `field`: alternativly a pre defined radiation field
- `water_attenuation_coefficient`: the radiation attenuation coefficient of the water
- `water_heat_capacity`: the specific heat capacity of the water
- `water_density`: density of the water

Example
=======

```jldoctest
julia> using Walrus: BodyHeating

julia> using Oceananigans

julia> grid = RectilinearGrid(size = (128, 128, 128), extent = (1000, 1000, 1000))
128×128×128 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 3×3×3 halo
├── Periodic x ∈ [0.0, 1000.0)  regularly spaced with Δx=7.8125
├── Periodic y ∈ [0.0, 1000.0)  regularly spaced with Δy=7.8125
└── Bounded  z ∈ [-1000.0, 0.0] regularly spaced with Δz=7.8125
julia> body_heating = BodyHeating(; surface_flux = (x, y, t) -> 100, grid)
┌ Warning: This radiative heating model is untested
└ @ Walrus.RadiativeTransfer ~/Documents/Projects/Walrus.jl/src/radiative_transfer/body_heating.jl:88
Single band light attenuation and body heating model with: 
- Attenuation coefficient: 1.8 (1 / m)
- Heat capacity: 3991.0 (J / K / kg)

julia> model = NonhydrostaticModel(; grid, biogeochemistry = body_heating, tracers = :T)
NonhydrostaticModel{CPU, RectilinearGrid}(time = 0 seconds, iteration = 0)
├── grid: 128×128×128 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 3×3×3 halo
├── timestepper: QuasiAdamsBashforth2TimeStepper
├── tracers: T
├── closure: Nothing
├── buoyancy: Nothing
└── coriolis: Nothing
```
"""
function BodyHeating(; surface_flux,
                      grid = nothing,
                      field = nothing,
                      water_attenuation_coefficient = 1.8,
                      water_heat_capacity = 3991.0, # J K⁻¹ kg⁻¹
                      water_density = 1026.0) # kg m⁻³

   @warn "This radiative heating model is untested"

   isnothing(grid) && isnothing(field) && error("You must specify either the radiation field or grid")

   if isnothing(field)
       field = RadiationField(grid, surface_flux)
   end

   return BodyHeating(water_attenuation_coefficient,
                      water_heat_capacity,
                      water_density,
                      surface_flux,
                      field)
end

summary(::BodyHeating) = string("Single band light attenuation and body heating model")
show(io::IO, body_heating::BodyHeating) = println(io, string(summary(body_heating), " with: \n",
                                                             "- Attenuation coefficient: ", body_heating.water_attenuation_coefficient, " (1 / m)\n",
                                                             "- Heat capacity: ", body_heating.water_heat_capacity, " (J / K / kg)"))

@kernel function update_body_heat!(body_heating, grid, t)
   i, j = @index(Global, NTuple)

   heating = body_heating.field

   α = body_heating.water_attenuation_coefficient

   x, y, _ = node(i, j, 1, grid, Center(), Center(), Center())
   
   surface_heating = body_heating.surface_flux(x, y, t)

   zᶠ = znodes(grid, Center(), Center(), Face())

   # the rest of the points
   @unroll for k in grid.Nz:-1:1
       @inbounds begin
           heating[i, j, k] = surface_heating * (exp(- α * abs(zᶠ[k + 1])) - exp(- α * abs(zᶠ[k])))
       end
   end
end

function update_biogeochemical_state!(model, body_heating::BodyHeating)
   arch = architecture(model.grid)

   launch!(arch, model.grid, :xy, update_body_heat!, body_heating, model.grid, model.clock.time)

   fill_halo_regions!(body_heating.field, model.clock, fields(model))
end

@kernel function apply_body_heating!(body_heating, grid, ∂ₜT)
   i, j, k = @index(Global, NTuple)

   heating = body_heating.field

   α = body_heating.water_attenuation_coefficient
   cᵖ = body_heating.water_heat_capacity
   ρₒ = body_heating.water_density

   @inbounds ∂ₜT[i, j, k] += α * heating[i, j, k] / (cᵖ * ρₒ)
end

function update_tendencies!(bgc, body_heating::BodyHeating, model)
   arch = architecture(model.grid)

   ∂ₜT = model.timestepper.Gⁿ.T

   launch!(arch, model.grid, :xyz, apply_body_heating!, body_heating, model.grid, ∂ₜT)
end