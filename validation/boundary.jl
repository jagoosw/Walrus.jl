using Walrus.SurfaceFluxModel: OceanAtmosphereBoundary, OceanAtmosphereBoundaryConditions, PrescribedAtmosphericState
using Oceananigans.Fields: ConstantField

atmosphere = PrescribedAtmosphericState(wind_speed = 2.0, wind_direction = 0.0, temperature = 15.0)

boundary_conditions = OceanAtmosphereBoundaryConditions(atmosphere)

boundary = boundary_conditions.u.condition.func

u = ConstantField(0)
v = ConstantField(0)
T = ConstantField(14.0)

grid = (; Nz = 1)

@info boundary(1, 1, grid, nothing, (; u, v, T), Val(:u))

@info boundary(1, 1, grid, nothing, (; u, v, T), Val(:v))

@info boundary(1, 1, grid, nothing, (; u, v, T), Val(:T))
