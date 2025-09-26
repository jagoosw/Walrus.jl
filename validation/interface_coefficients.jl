using CairoMakie

using Walrus.SurfaceFluxModel: SimilarityTheoryInterface, SmoothAndCharnock

interface = SimilarityTheoryInterface(; roughness_length = SmoothAndCharnock(air_kinematic_viscosity=14e-6))

Us = Float64[2, 5, 10, 15, 20, 25]#[0:0.2:25;]
ΔTs = Float64[-20:5:-5..., -4:5..., 10:5:20...]#Float64[[-20:0;]..., [1:5;]..., [10:5:40;]...]#[-0.1:0.001:0.1;]#
T₀ = Float64(12)

Cₘ = zeros(length(Us), length(ΔTs))
Cₕ = zeros(length(Us), length(ΔTs))

τ = zeros(length(Us), length(ΔTs))
Q = zeros(length(Us), length(ΔTs))

for (i, U) in enumerate(Us), (j, ΔT) in enumerate(ΔTs)
    @info U, ΔT
    Cₘ[i, j], Cₕ[i, j] = interface(U, T₀, T₀ + ΔT + 0.1, 0)

    τ[i, j] = Cₘ[i, j] * U^2
    Q[i, j] = Cₕ[i, j] * U * ΔT
end