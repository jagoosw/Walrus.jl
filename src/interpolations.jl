module Interpolations

export SimpleInterpolation

using Adapt: adapt

using Oceananigans.Architectures: on_architecture, CPU

import Adapt: adapt_structure

struct SimpleInterpolation{R, V, M} <: Function
               range :: R
              values :: V
  boundary_condition :: M
end

adapt_structure(to, itp::SimpleInterpolation) = SimpleInterpolation(adapt(to, itp.range),
                                                                    adapt(to, itp.values))

function coordinate_to_range(coordinate)
    x₀ = minimum(coordinate)
    x₁ = maximum(coordinate)

    (coordinate[2] - coordinate[1] ≈ coordinate[end] - coordinate[end - 1]) || throw(ArgumentError("Interpolation coordinate must be regularly spaced"))

    Δx = coordinate[2] - coordinate[1]
    
    return (; x₀, x₁, Δx)
end

function SimpleInterpolation(coordinate::Array, values; arch = CPU(), boundary_condition = Cyclic())
    range = coordinate_to_range(coordinate)

    return SimpleInterpolation(range, on_architecture(arch, values), boundary_condition)
end

function SimpleInterpolation(coordinates::Tuple, values; arch = CPU(), boundary_conditions = Cyclic())
    if !(boundary_conditions isa Tuple)
        boundary_conditions = tuple([boundary_conditions for n in 1:length(coordinates)]...)
    end

    length(coordinates) == length(boundary_conditions) || 
        throw(ArgumentError("You must either provide one boundary condition to be used in all directions, or one for each direction"))

    ranges = tuple([coordinate_to_range(coord) for coord in coordinates])

    return SimpleInterpolation(ranges, on_architecture(arch, values), boundary_condition)
end

(itp::SimpleInterpolation)(X...) = itp(X..., itp.values, itp.range, itp.boundary_condition)

struct Limited end
struct Cyclic end

@inline function (::Limited)(x, x₀, Δx, N)
    # this will (hopeflly) result in an out of bounds error 
    # but the whole point of this excersise was to make it kernel safe
    # so I'm not erroring here
    n₀ = floor(Int, (x - x₀) / Δx)

    return x, n₀ + 1, n₀ + 2
end

@inline function (::Cyclic)(x, x₀, Δx, N)
    n₀ = floor(Int, (x - x₀) / Δx)

    n₁ = mod(n₀, N) + 1

    n₂ = ifelse(n₁ == N, 1, n₁ + 1)

    x = mod(x - x₀, N * Δx) + x₀

    return x, n₁, n₂
end

# the terminal case
function (itp::SimpleInterpolation)(x::Number, values, range, boundary_condition)
    x₀ = range.x₀
    x₁ = range.x₁
    Δx = range.Δx
    N  = length(values)

    x, n₁, n₂ = boundary_condition(x, x₀, Δx, N)

    x₁ = x₀ + Δx * (n₁ - 1)

    y₁ = @inbounds values[n₁]
    y₂ = @inbounds values[n₂]

    return y₁ + (x - x₁) * (y₂ - y₁) / Δx
end

end # module