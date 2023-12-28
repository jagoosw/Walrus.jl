module Interpolations

export SimpleInterpolation

using Adapt: adapt

import Adapt: adapt_structure

struct SimpleInterpolation{R, V}
     range :: R
    values :: V
end

adapt_structure(to, itp::SimpleInterpolation) = SimpleInterpolation(adapt(to, itp.range),
                                                                    adapt(to, itp.values))


function (itp::SimpleInterpolation)(x)
    n₁ = floor(Int, (x - itp.range.x₀) / itp.range.Δx)

    x₁ = itp.range.x₀ + itp.range.Δx * n₁
    x₂ = x₁ + itp.range.Δx

    y₁ = @inbounds itp.values[n₁ + 1]
    y₂ = @inbounds itp.values[n₁ + 2]

    return y₁ + (x - x₁) * (y₂ - y₁) / (x₂ - x₁)
end

end # module