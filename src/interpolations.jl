module Interpolations

export SimpleInterpolation

using Adapt: adapt

using Oceananigans.Architectures: on_architecture, CPU

import Adapt: adapt_structure

struct SimpleInterpolation{R, V}
     range :: R
    values :: V
end

adapt_structure(to, itp::SimpleInterpolation) = SimpleInterpolation(adapt(to, itp.range),
                                                                    adapt(to, itp.values))

function SimpleInterpolation(range::Array, values; arch = CPU())
    x₀ = minimum(range)

    (range[2] - range[1] ≈ range[end] - range[end - 1]) || throw(ArgumentError("Interpolation range must be regularly spaced"))

    Δx = range[2] - range[1]

    return SimpleInterpolation((; x₀, Δx), on_architecture(arch, values))
end


function (itp::SimpleInterpolation)(x)
    n₁ = floor(Int, (x - itp.range.x₀) / itp.range.Δx)

    x₁ = itp.range.x₀ + itp.range.Δx * n₁
    x₂ = x₁ + itp.range.Δx

    y₁ = @inbounds itp.values[n₁ + 1]
    y₂ = @inbounds itp.values[n₁ + 2]

    return y₁ + (x - x₁) * (y₂ - y₁) / (x₂ - x₁)
end

end # module