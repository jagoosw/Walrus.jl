"""
    HomogeneousBodyHeating

A model for single band light attenuation which heats the water in the form:
```math
I(x, y, z) = I_0(x, y) * \\exp\\left(-\\alpha z\\right),
```
where ``I`` is the radiation intensity and ``\\alpha`` is the attenuation coefficient. This heats the water
like ``\\frac{\\partial T(x, y, z)}{\\partial t} = \\frac{I(x, y, z)A}{c^p\\rho}`` where ``A`` is the area of the cell, ``c^p`` is the specific heat capacity
and ``\\rho`` is the water density. 
"""
struct HomogeneousBodyHeating{FT, S} <: Function
 water_attenuation_coefficient :: FT
           water_heat_capacity :: FT
                 water_density :: FT

                  surface_flux :: S
end

adapt_structure(to, bh::HomogeneousBodyHeating) = 
    HomogeneousBodyHeating(bh.water_attenuation_coefficient,
                           bh.water_heat_capacity,
                           bh.water_density,
                           
                           adapt(to, bh.surface_flux))

"""
    HomogeneousBodyHeating(; surface_flux,
                             water_attenuation_coefficient = 1.8,
                             water_heat_capacity = 3991.0, # J K⁻¹ kg⁻¹
                             water_density = 1026.0) # kg m⁻³

Creates a model in which a `surface_flux` (W / m²) is attenuated by and heats the water.
This interacts with Oceananigans as a body forcing.

Keyword Arguments
=================

- `surface_flux` (required): a function returning the surface radiaiton flux in the form `surface_flux(x, y, t)` or single value
- `water_attenuation_coefficient`: the radiation attenuation coefficient of the water
- `water_heat_capacity`: the specific heat capacity of the water
- `water_density`: density of the water

Example
=======

```jldoctest; filter = r".*@ Walrus.RadiativeTransfer.*"
julia> using Walrus: HomogeneousBodyHeating

julia> using Oceananigans

julia> grid = RectilinearGrid(size = (128, 128, 128), extent = (1000, 1000, 1000))
128×128×128 RectilinearGrid{Float64, Oceananigans.Grids.Periodic, Oceananigans.Grids.Periodic, Oceananigans.Grids.Bounded} on Oceananigans.Architectures.CPU with 3×3×3 halo
├── Periodic x ∈ [0.0, 1000.0)  regularly spaced with Δx=7.8125
├── Periodic y ∈ [0.0, 1000.0)  regularly spaced with Δy=7.8125
└── Bounded  z ∈ [-1000.0, 0.0] regularly spaced with Δz=7.8125
julia> body_heating = HomogeneousBodyHeating(; surface_flux = (x, y, t) -> 100)
(::HomogeneousBodyHeating{Float64, Walrus.ContinuousSurfaceFunction{var"#1#2"}}) (generic function with 1 method)

julia> model = NonhydrostaticModel(; grid, forcing = (; T = Forcing(body_heating, discrete_form=true)), tracers = :T)
NonhydrostaticModel{CPU, RectilinearGrid}(time = 0 seconds, iteration = 0)
├── grid: 128×128×128 RectilinearGrid{Float64, Oceananigans.Grids.Periodic, Oceananigans.Grids.Periodic, Oceananigans.Grids.Bounded} on Oceananigans.Architectures.CPU with 3×3×3 halo
├── timestepper: RungeKutta3TimeStepper
├── advection scheme: Centered(order=2)
├── tracers: T
├── closure: Nothing
├── buoyancy: Nothing
└── coriolis: Nothing
```
"""
function HomogeneousBodyHeating(; surface_flux,
                                  water_attenuation_coefficient = 1.8,
                                  water_heat_capacity = 3991.0, # J K⁻¹ kg⁻¹
                                  water_density = 1026.0) # kg m⁻³

    surface_flux = normalise_surface_function(surface_flux)

    return HomogeneousBodyHeating(water_attenuation_coefficient,
                                  water_heat_capacity,
                                  water_density,
                                  surface_flux)
end

@inline function (heating::HomogeneousBodyHeating)(i, j, k, grid, clock, model_fields)
    ρ = heating.water_density
    cₚ = heating.water_heat_capacity
    α = heating.water_attenuation_coefficient

    x, y, _ = node(i, j, k, grid, Center(), Center(), Center())

    zᶠ = znode(i, j, k, grid, Center(), Center(), Face())

    zᶠ⁺ = znode(i, j, k + 1, grid, Center(), Center(), Face())

    return get_value(heating.surface_flux, i, j, grid, clock) * (exp(- α * abs(zᶠ⁺)) - exp(- α * abs(zᶠ))) / (ρ * cₚ) / (zᶠ⁺ - zᶠ)
end


summary(::HomogeneousBodyHeating) = string("Single band light attenuation and body heating model")
show(io::IO, body_heating::HomogeneousBodyHeating) = println(io, string(summary(body_heating), " with: \n",
                                                             " Attenuation coefficient: ", body_heating.water_attenuation_coefficient, " (1 / m)\n",
                                                             " Heat capacity: ", body_heating.water_heat_capacity, " (J / K / kg)"))
