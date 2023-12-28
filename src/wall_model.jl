"""
WallStressModel

A wall stress model for LES simulation with default parameters similar
to that proposed in [SCHUMANN1975376](@citet), [HARTEL1996283](@citet),
[Piomelli1989](@citet), and [taylor2007](@citet).
"""    
module WallStressModel

export WallStress, WallStressBoundaryConditions

using Roots

using Adapt: adapt

using Oceananigans.BoundaryConditions: FluxBoundaryCondition
using Oceananigans.Fields: Center, Face
using Oceananigans.Grids: znode

using Oceananigans.Architectures: arch_array, CPU, architecture

using Walrus.Interpolations: SimpleInterpolation

import Adapt: adapt_structure
import Base: summary, show

struct WallStress{FT, U} <: Function
    von_Karman_constant :: FT
    kinematic_viscosity :: FT
                      B :: FT

      friction_velocities :: U
end

adapt_structure(to, ws::WallStress) = WallStress(ws.von_Karman_constant, ws.kinematic_viscosity,
                                                 ws.B, adapt(to, ws.friction_velocities))

"""
    WallStress(; von_Karman_constant = 0.4,
                 kinematic_viscosity = 1e-6,
                 B = 5.2,
                 precomputed_friction_velocities = false,
                 precompute_speeds = [0:25/100000:25;],
                 grid = nothing)

Returns a wall stress model for LES simulation with default parameters similar
to that proposed in [SCHUMANN1975376](@citet), [HARTEL1996283](@citet),
[Piomelli1989](@citet), and [taylor2007](@citet).

Friction velocities will be precomputed at `precompute_speeds` if 
`precomputed_friction_velocities` is true and `grid` is provided.

Keyword Arguments
=================

- `von_Karman_constant`: the von Karman wall stress constant
- `kinematic_viscosity`: kinematic viscosity of the water above the wall
- `B`: wall stress constant 
- `precomputed_friction_velocities`: precompute friction velocities?
- `precompute_speeds`: bottom water speeds to precompute friction velocities for, 
  this should encompas the range of speeds possible in your simulation
- `grid`: the grid to precompute the friction velocities for

Example
=======

```jldoctest
julia> using Walrus: WallStress

julia> using Oceananigans

julia> wall_stress = WallStress()
(::WallStress{Float64, Nothing}) (generic function with 1 method)
julia> boundary_conditions = (u = FieldBoundaryConditions(bottom = FluxBoundaryCondition(wall_stress, discrete_form = true, parameters = Val(:x))),
                              v = FieldBoundaryConditions(bottom = FluxBoundaryCondition(wall_stress, discrete_form = true, parameters = Val(:y))))
(u = Oceananigans.FieldBoundaryConditions, with boundary conditions
├── west: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── east: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── south: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── north: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── bottom: FluxBoundaryCondition: DiscreteBoundaryFunction (::WallStress{Float64, Nothing}) with parameters Val{:x}
├── top: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
└── immersed: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing), v = Oceananigans.FieldBoundaryConditions, with boundary conditions
├── west: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── east: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── south: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── north: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── bottom: FluxBoundaryCondition: DiscreteBoundaryFunction (::WallStress{Float64, Nothing}) with parameters Val{:y}
├── top: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
└── immersed: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing))

```
"""
function WallStress(; von_Karman_constant::FT = 0.4,
                      kinematic_viscosity::FT = 1.15e-6,
                      B::FT = 5.2,
                      precomputed_friction_velocities = false,
                      precompute_speeds = [0:25/100000:25;],
                      grid = nothing,
                      arch = isnothing(grid) ? CPU() : architecture(grid)) where FT
    
    if precomputed_friction_velocities
        tmp = WallStress(von_Karman_constant, kinematic_viscosity, B, nothing)

        ν = kinematic_viscosity
        κ = von_Karman_constant

        z₁ = znode(1, grid, Center()) - znode(1, grid, Face())

        n = 100000
        velocities = zeros(length(precompute_speeds))

        for (n, speed) in enumerate(precompute_speeds)
            params = (; z₁, κ, ν, B)
            velocities[n] = find_friction_velocity(tmp, speed, params)
        end

        friction_velocities = SimpleInterpolation(precompute_speeds, velocities; arch)
    else
        friction_velocities = nothing
    end

    return WallStress(von_Karman_constant, kinematic_viscosity, B, friction_velocities)
end

@inline stress_velocity(uₜ, params) = log(params.z₁ * uₜ / (params.ν + eps(0.0))) / (params.κ + eps(0.0)) + params.B - params.U₁ / (uₜ + eps(0.0))

@inline function find_friction_velocity(::WallStress{<:Any, Nothing}, U₁, params)
    uₜ = 0
    
    U₁ == 0 || (uₜ = find_zero(stress_velocity, (0., Inf), Bisection(); p = merge(params, (; U₁)), maxiters = 10^5))

    return uₜ
end

@inline find_friction_velocity(ws::WallStress, U₁, params) = 
    ws.friction_velocities(U₁)

@inline function (wall_stress::WallStress)(i, j, grid, clock, model_fields, ::Val{direction}) where direction
    ν = wall_stress.kinematic_viscosity
    κ = wall_stress.von_Karman_constant
    B = wall_stress.B

    @inbounds begin
        u = model_fields.u[i, j, 1]
        v = model_fields.v[i, j, 1]

        z₁ = znode(i, j, 1, grid, Center(), Center(), Center()) - znode(i, j, 1, grid, Center(), Center(), Face())
    end

    U₁ = sqrt(u ^ 2 + v ^ 2)

    uₜ = find_friction_velocity(wall_stress, U₁, (; z₁, κ, ν, B))

    return ifelse(direction == :x, -u / (U₁ + eps(0.0)) * uₜ ^ 2, -v / (U₁ + eps(0.0)) * uₜ ^ 2)
end


"""
    WallStressBoundaryConditions(; von_Karman_constant = 0.4,
                                   kinematic_viscosity = 1e-6,
                                   B = 5.2,
                                   precompute_speeds = [0:25/100000:25;],
                                   grid = nothing)

Convenience constructor to setup `WallStress` boundary conditions.

Keyword Arguments
=================

- `von_Karman_constant`: the von Karman wall stress constant
- `kinematic_viscosity`: kinematic viscosity of the water above the wall
- `B`: wall stress constant 
- `precomputed_friction_velocities`: precompute friction velocities?
- `precompute_speeds`: bottom water speeds to precompute friction velocities for, 
  this should encompas the range of speeds possible in your simulation
- `grid`: the grid to precompute the friction velocities for

Example
=======

```jldoctest
julia> using Walrus: WallStressBoundaryConditions

julia> using Oceananigans

julia> stress_boundary_conditions = WallStressBoundaryConditions()
(u = FluxBoundaryCondition: DiscreteBoundaryFunction (::WallStress{Float64, Nothing}) with parameters Val{:x}, v = FluxBoundaryCondition: DiscreteBoundaryFunction (::WallStress{Float64, Nothing}) with parameters Val{:y})
julia> boundary_conditions = (u = FieldBoundaryConditions(bottom = stress_boundary_conditions.u),
                              v = FieldBoundaryConditions(bottom = stress_boundary_conditions.v))
(u = Oceananigans.FieldBoundaryConditions, with boundary conditions
├── west: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── east: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── south: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── north: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── bottom: FluxBoundaryCondition: DiscreteBoundaryFunction (::WallStress{Float64, Nothing}) with parameters Val{:x}
├── top: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
└── immersed: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing), v = Oceananigans.FieldBoundaryConditions, with boundary conditions
├── west: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── east: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── south: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── north: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
├── bottom: FluxBoundaryCondition: DiscreteBoundaryFunction (::WallStress{Float64, Nothing}) with parameters Val{:y}
├── top: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)
└── immersed: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing))
```
"""
function WallStressBoundaryConditions(; von_Karman_constant::FT = 0.4,
                                        kinematic_viscosity::FT = 1.15e-6,
                                        B::FT = 5.2,
                                        precomputed_friction_velocities = false,
                                        precompute_speeds = [0:25/100000:25;],
                                        grid = nothing) where FT

    wall_stress_instance = WallStress(; von_Karman_constant,
                                        kinematic_viscosity,
                                        B,
                                        precomputed_friction_velocities,
                                        precompute_speeds,
                                        grid)

    u = FluxBoundaryCondition(wall_stress_instance, discrete_form = true, parameters = Val(:x))

    v = FluxBoundaryCondition(wall_stress_instance, discrete_form = true, parameters = Val(:y))

    return (; u, v)
end

summary(::WallStress) = string("Wall stress model")
show(io::IO, wall_stress::WallStress) = println(io, "Wall stress model with ν = $(wall_stress.kinematic_viscosity), κ = $(wall_stress.von_Karman_constant), and B = $(wall_stress.B)")

end # module