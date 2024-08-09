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

@inline (::Linear)(x, x₀, x₁) = x
@inline (::Cyclic)(x, x₀, x₁) = mod(x - x₀, x₁ - x₀) + x₀

function (itp::SimpleInterpolation)(x)
    x = itp.mode(x, itp.range.x₀, itp.range.x₁)

    n₁ = floor(Int, (x - itp.range.x₀) / itp.range.Δx)

    x₁ = itp.range.x₀ + itp.range.Δx * n₁
    x₂ = x₁ + itp.range.Δx

    y₁ = @inbounds itp.values[n₁ + 1]
    y₂ = @inbounds itp.values[n₁ + 2]

    return y₁ + (x - x₁) * (y₂ - y₁) / (x₂ - x₁)
end

end # module