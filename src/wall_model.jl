"""
WallStressModel

A wall stress model for LES simulation with default parameters similar
to that proposed in [SCHUMANN1975376](@citet), [HARTEL1996283](@citet),
[Piomelli1989](@citet), and [taylor2007](@citet).
"""    
module WallStressModel

export WallStress, WallStressBoundaryConditions

using Roots

using Oceananigans.BoundaryConditions: FluxBoundaryCondition
using Oceananigans.Fields: Center, Face
using Oceananigans.Grids: znode

"""
    WallStress(; von_Karman_constant = 0.4,
                 kinematic_viscosity = 1e-6,
                 B = 5.2)

Returns a wall stress model for LES simulation with default parameters similar
to that proposed in [SCHUMANN1975376](@citet), [HARTEL1996283](@citet),
[Piomelli1989](@citet), and [taylor2007](@citet).


Keyword Arguments
=================

- `von_Karman_constant`: the von Karman wall stress constant
- `kinematic_viscosity`: kinematic viscosity of the water above the wall
- `B`: wall stress constant 

Example
=======

```jldoctest
julia> using Walrus: WallStressModel

julia> using Oceananigans

julia> wall_stress = WallStress()

julia> boundary_conditions = (u = FieldBoundaryConditions(bottom = FluxBoundaryCondition(wall_stress, discrete_form = true, parameters = :x)),
                              v = FieldBoundaryConditions(bottom = FluxBoundaryCondition(wall_stress, discrete_form = true, parameters = :y)))

```
"""

@kwdef struct WallStress{FT} <: Function
    von_Karman_constant :: FT = 0.4
    kinematic_viscosity :: FT = 1e-6
                      B :: FT = 5.2
end

@inline stress_velocity(uₜ, params) = log(params.z₁ * uₜ / (params.ν + eps(0.0))) / (params.κ + eps(0.0)) + params.B - params.U₁ / (uₜ + eps(0.0))

@inline function (wall_stress::WallStress)(i, j, grid, clock, model_fields, direction)
    ν = wall_stress.kinematic_viscosity
    κ = wall_stress.von_Karman_constant
    B = wall_stress.B

    @inbounds begin
        u = abs(model_fields.u[i, j, 1])
        v = abs(model_fields.v[i, j, 1])

        z₁ = znode(i, j, 1, grid, Center(), Center(), Center()) - znode(i, j, 1, grid, Center(), Center(), Face())
    end

    U₁ = sqrt(u ^ 2 + v ^ 2)

    uₜ = find_zero(stress_velocity, (0., Inf), Bisection(); p = (; U₁, z₁, κ, ν, B), maxiters = 10^5)

    return ifelse(direction == :x, -u / (U₁ + eps(0.0)) * uₜ ^ 2, -v / (U₁ + eps(0.0)) * uₜ ^ 2)
end


"""
    WallStressBoundaryConditions(; von_Karman_constant = 0.4,
                                   kinematic_viscosity = 1e-6,
                                   B = 5.2)

Convenience constructor to setup `WallStress` boundary conditions.

Keyword Arguments
=================

- `von_Karman_constant`: the von Karman wall stress constant
- `kinematic_viscosity`: kinematic viscosity of the water above the wall
- `B`: wall stress constant 

Example
=======

```jldoctest
julia> using Walrus: WallStressModel

julia> using Oceananigans

julia> stress_boundary_conditions = WallStressBoundaryConditions()

julia> boundary_conditions = (u = FieldBoundaryConditions(bottom = stress_boundary_conditions.u),
                              v = FieldBoundaryConditions(bottom = stress_boundary_conditions.v))

```
"""
function WallStressBoundaryConditions(; von_Karman_constant = 0.4,
                                        kinematic_viscosity = 1e-6,
                                        B = 5.2)

    wall_stress_instance = WallStress(von_Karman_constant,
                                      kinematic_viscosity,
                                      B)

    u = FluxBoundaryCondition(wall_stress_instance, discrete_form = true, parameters = :x)

    v = FluxBoundaryCondition(wall_stress_instance, discrete_form = true, parameters = :y)

    return (; u, v)
end

end # module