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

@inline function (::Linear)(x, x₀, x₁, Δx, N)
    n₀ = floor(Int, (x - x₀) / Δx)

    return n₀ + 1, n₀ + 2
end

@inline function (::Cyclic)(x, x₀, x₁, Δx, N)
    x = mod(x - x₀, x₁ - x₀ + Δx) + x₀

    n₀ = floor(Int, (x - x₀) / Δx)

    n₁, n₂ = n₀ + 1, n₀ + 2

    n₁, n₂ = ifelse(x > x₁, (N, 1), (n₁, n₂))

    return n₁, n₂
end

function (itp::SimpleInterpolation)(x)
    x = itp.mode(x, itp.range.x₀, itp.range.x₁)

    n₁, n₂ = floor(Int, (x - itp.range.x₀) / itp.range.Δx)

    x₁ = itp.range.x₀ + itp.range.Δx * (n₁ - 1)

    y₁ = @inbounds itp.values[n₁]
    y₂ = @inbounds itp.values[n₂]

    return y₁ + (x - x₁) * (y₂ - y₁) / itp.range.Δx
end

end # module