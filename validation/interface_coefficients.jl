using CairoMakie, Oceananigans
using Oceananigans.Fields: ConstantField, ZeroField

using Walrus.SurfaceFluxModel: SimilarityTheoryInterface, SmoothAndCharnock, update_interface!, PrescribedAtmosphericState

grid = RectilinearGrid(size = (1, 1, 1), extent = (1, 1, 1))

U = Ref(0.0)
T = Ref(0.0)

atmosphere = PrescribedAtmosphericState(wind_speed = (args...)->U[], wind_direction = 90.0, temperature = (args...) -> T[])

interface = SimilarityTheoryInterface(grid; roughness_length = SmoothAndCharnock(air_kinematic_viscosity=14e-6))

Us = Float64[2, 5, 10, 15, 20, 25]#[0:0.2:25;]
ΔTs = Float64[-20:5:-5..., -4:5..., 10:5:20...]#Float64[[-20:0;]..., [1:5;]..., [10:5:40;]...]#[-0.1:0.001:0.1;]#
T₀ = Float64(12)

Cₘ = zeros(length(Us), length(ΔTs))
Cₕ = zeros(length(Us), length(ΔTs))

τ = zeros(length(Us), length(ΔTs))
Q = zeros(length(Us), length(ΔTs))

import Oceananigans: fields

fields(model::NamedTuple) = model.fields

for (i, U1) in enumerate(Us), (j, ΔT) in enumerate(ΔTs)
    @info U1, ΔT

    U[] = U1
    T[] = T₀ + ΔT

    update_interface!(interface, (; clock = (; time = 0.0), grid, fields = (; T = ConstantField(T₀), u = ZeroField(), v = ZeroField())), atmosphere)

    Cₘ[i, j], Cₕ[i, j] = interface.drag_coefficient[1, 1, 1], interface.heat_exchange_coefficient[1, 1, 1]

    τ[i, j] = Cₘ[i, j] * U1^2
    Q[i, j] = Cₕ[i, j] * U1 * ΔT
end