module Interpolations

export SimpleInterpolation

using Adapt: adapt

using Oceananigans.Architectures: on_architecture, CPU

import Adapt: adapt_structure

struct SimpleInterpolation{R, V, C} <: Function
     range :: R
    values :: V
   closure :: C
end

adapt_structure(to, itp::SimpleInterpolation) = SimpleInterpolation(adapt(to, itp.range),
                                                                    adapt(to, itp.values),
                                                                    adapt(to, itp.closure))

function SimpleInterpolation(range::Array, values; closure = Cyclic(length(range)), arch = CPU())
    x₀ = minimum(range)

    (range[2] - range[1] ≈ range[end] - range[end - 1]) || throw(ArgumentError("Interpolation range must be regularly spaced"))

    Δx = range[2] - range[1]

    return SimpleInterpolation((; x₀, Δx), on_architecture(arch, values), closure)
end

function (itp::SimpleInterpolation)(x)
    n₁, n₂, x₁ = get_indices(itp.closure, x, itp.range)

    y₁ = @inbounds itp.values[n₁]
    y₂ = @inbounds itp.values[n₂]

    return y₁ + (x - x₁) * (y₂ - y₁) / itp.range.Δx
end

@inline (itp::SimpleInterpolation)(x, y, t) = itp(t) 

struct Cyclic{IT} 
    n :: IT

    function Cyclic(n::IT) where IT
        @warn "You are using a cyclic closure, make sure that the final value in your range corresponds to the start value - Δx"

        return new{IT}(n)
    end
end

@inline function get_indices(c::Cyclic, x, range)
    N = mod((x - range.x₀) / range.Δx, c.n)

    n₁, n₂ = ifelse(N+1 >= c.n, (c.n, 1), (floor(Int, N) + 1, floor(Int, N) + 2))

    x₁ = range.x₀ + range.Δx * (floor(Int, (x - range.x₀) / range.Δx))

    return n₁, n₂, x₁
end

end # module
