module Interpolations

export SimpleInterpolation

using Adapt: adapt

using Oceananigans.Architectures: on_architecture, CPU

import Adapt: adapt_structure

struct Linear end
struct Cyclic end

struct SimpleInterpolation{R, V, M} <: Function
     range :: R
    values :: V
      mode :: M
end

adapt_structure(to, itp::SimpleInterpolation) = SimpleInterpolation(adapt(to, itp.range),
                                                                    adapt(to, itp.values))

function SimpleInterpolation(range::Array, values; arch = CPU(), mode = Cyclic())
    x₀ = minimum(range)
    x₁ = maximum(range)

    (range[2] - range[1] ≈ range[end] - range[end - 1]) || throw(ArgumentError("Interpolation range must be regularly spaced"))

    Δx = range[2] - range[1]

    return SimpleInterpolation((; x₀, x₁, Δx), on_architecture(arch, values), mode)
end

@inline function (::Linear)(x, x₀, Δx, N)
    n₀ = floor(Int, (x - x₀) / Δx)

    return x, n₀ + 1, n₀ + 2
end

@inline function (::Cyclic)(x, x₀, Δx, N)
    n₀ = floor(Int, (x - x₀) / Δx)

    n₁ = mod(n₀, N - 1) + 1

    n₂ = ifelse(n₁ == N, 1, n₁ + 1)

    x = x₀ + Δx * (n₁ - 1)

    return x, n₁, n₂
end

function (itp::SimpleInterpolation)(x)
    x₀ = itp.range.x₀
    x₁ = itp.range.x₁
    Δx = itp.range.Δx
    N  = length(itp.values)

    x, n₁, n₂ = itp.mode(x, x₀, Δx, N)

    x₁ = x₀ + Δx * (n₁ - 1)

    y₁ = @inbounds itp.values[n₁]
    y₂ = @inbounds itp.values[n₂]

    return y₁ + (x - x₁) * (y₂ - y₁) / Δx
end

end # module