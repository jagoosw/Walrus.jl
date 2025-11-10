module Interpolations

export SimpleInterpolation

using Adapt: adapt

using Oceananigans.Architectures: on_architecture, CPU

import Adapt: adapt_structure

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

struct Fixed{IT} 
    n :: IT
end

@inline function get_indices(c::Fixed, x, range)
    N = (x - range.x₀) / range.Δx

    n₁, n₂ = floor(Int, N) + 1, floor(Int, N) + 2

    n₁ = ifelse(n₁ > c.n, c.n, n₁)
    n₂ = ifelse(n₂ > c.n, c.n, n₂)

    n₁ = ifelse(n₁ < 1, 1, n₁)
    n₂ = ifelse(n₂ < 1, 1, n₂)

    x₁ = range.x₀ + range.Δx * (floor(Int, (x - range.x₀) / range.Δx))

    return n₁, n₂, x₁
end

#####
##### New interpolation
#####

struct MultiDimensionalInterpolation{N, D, C} <: Function
           data :: D
    coordinates :: C
end

function adapt_structure(to, itp::MultiDimensionalInterpolation{N}) where N 
    data = adapt(to, itp.data)
    coordinates = adapt(to, itp.coordinates)

    return MultiDimensionalInterpolation{N, typeof(data), typeof(coordinates)}(data, coordinates)
end

function MultiDimensionalInterpolation(ranges, values, arch = CPU())
    coordinates = tuple(map(
        n -> (Δx = diff(ranges[n])[1], x₀ = minimum(ranges[n]), closure = Fixed(length(ranges[n]))),
        1:length(ranges)
    )...)

    data = on_architecture(arch, values)

    return MultiDimensionalInterpolation{length(ranges), typeof(data), typeof(coordinates)}(data, coordinates)
end

function MultiDimensionalInterpolation(range::Union{Array, StepRange, StepRangeLen}, values; closure = Cyclic(length(range)), arch = CPU())
    x₀ = minimum(range)

    (range[2] - range[1] ≈ range[end] - range[end - 1]) || throw(ArgumentError("Interpolation range must be regularly spaced"))

    Δx = range[2] - range[1]
    coordinates = ((; x₀, Δx, closure), )

    data = on_architecture(arch, values)

    return MultiDimensionalInterpolation{1, typeof(data), typeof(coordinates)}(data, coordinates)
end

@inline function (itp::MultiDimensionalInterpolation{1})(x)
    crd = itp.coordinates[1]
    n₁, n₂, x₁ = get_indices(crd.closure, x, crd)

    y₁ = @inbounds itp.data[n₁]
    y₂ = @inbounds itp.data[n₂]

    return y₁ + (x - x₁) * (y₂ - y₁) / crd.Δx
end

# hmmmmmmm
@inline (itp::MultiDimensionalInterpolation{1})(x, y, t) = itp(t)

@inline function (itp::MultiDimensionalInterpolation{2})(x, y)
    crd = itp.coordinates

    @inbounds begin
        i₁, i₂, x₁ = get_indices(crd[1].closure, x, crd[1])
        j₁, j₂, y₁ = get_indices(crd[2].closure, y, crd[2])

        c₁₁ = itp.data[i₁, j₁]
        c₁₂ = itp.data[i₁, j₂]
        c₂₁ = itp.data[i₂, j₁]
        c₂₂ = itp.data[i₂, j₂]

        c₁ = c₁₁ + (y - y₁) * (c₁₂ - c₁₁) / crd[2].Δx
        c₂ = c₂₁ + (y - y₁) * (c₂₂ - c₂₁) / crd[2].Δx

        return c₁ + (x - x₁) * (c₂ - c₁) / crd[1].Δx
    end
end

@inline function (itp::MultiDimensionalInterpolation{3})(x, y, z)
    crd = itp.coordinates

    @inbounds begin
        i₁, i₂, x₁ = get_indices(crd[1].closure, x, crd[1])
        j₁, j₂, y₁ = get_indices(crd[2].closure, y, crd[2])
        k₁, k₂, z₁ = get_indices(crd[3].closure, z, crd[3])

        c₁₁₁ = itp.data[i₁, j₁, k₁]
        c₁₂₁ = itp.data[i₁, j₂, k₁]
        c₂₁₁ = itp.data[i₂, j₁, k₁]
        c₂₂₁ = itp.data[i₂, j₂, k₁]
        c₁₁₂ = itp.data[i₁, j₁, k₂]
        c₁₂₂ = itp.data[i₁, j₂, k₂]
        c₂₁₂ = itp.data[i₂, j₁, k₂]
        c₂₂₂ = itp.data[i₂, j₂, k₂]

        c₁₁ = c₁₁₁ + (z - z₁) * (c₁₁₂ - c₁₁₁) / crd[3].Δx
        c₂₁ = c₂₁₁ + (z - z₁) * (c₂₁₂ - c₂₁₁) / crd[3].Δx
        c₁₂ = c₁₂₁ + (z - z₁) * (c₁₂₂ - c₁₂₁) / crd[3].Δx
        c₂₂ = c₂₂₁ + (z - z₁) * (c₂₂₂ - c₂₂₁) / crd[3].Δx

        c₁ = c₁₁ + (y - y₁) * (c₁₂ - c₁₁) / crd[2].Δx
        c₂ = c₂₁ + (y - y₁) * (c₂₂ - c₂₁) / crd[2].Δx

        return c₁ + (x - x₁) * (c₂ - c₁) / crd[1].Δx
    end
end

@inline function (itp::MultiDimensionalInterpolation{4})(x, y, z, w)
    crd = itp.coordinates

    @inbounds begin
        i₁, i₂, x₁ = get_indices(crd[1].closure, x, crd[1])
        j₁, j₂, y₁ = get_indices(crd[2].closure, y, crd[2])
        k₁, k₂, z₁ = get_indices(crd[3].closure, z, crd[3])
        l₁, l₂, w₁ = get_indices(crd[4].closure, w, crd[4])

        c₁₁₁₁ = itp.data[i₁, j₁, k₁, l₁]
        c₁₂₁₁ = itp.data[i₁, j₂, k₁, l₁]
        c₂₁₁₁ = itp.data[i₂, j₁, k₁, l₁]
        c₂₂₁₁ = itp.data[i₂, j₂, k₁, l₁]
        c₁₁₂₁ = itp.data[i₁, j₁, k₂, l₁]
        c₁₂₂₁ = itp.data[i₁, j₂, k₂, l₁]
        c₂₁₂₁ = itp.data[i₂, j₁, k₂, l₁]
        c₂₂₂₁ = itp.data[i₂, j₂, k₂, l₁]
        c₁₁₁₂ = itp.data[i₁, j₁, k₁, l₂]
        c₁₂₁₂ = itp.data[i₁, j₂, k₁, l₂]
        c₂₁₁₂ = itp.data[i₂, j₁, k₁, l₂]
        c₂₂₁₂ = itp.data[i₂, j₂, k₁, l₂]
        c₁₁₂₂ = itp.data[i₁, j₁, k₂, l₂]
        c₁₂₂₂ = itp.data[i₁, j₂, k₂, l₂]
        c₂₁₂₂ = itp.data[i₂, j₁, k₂, l₂]
        c₂₂₂₂ = itp.data[i₂, j₂, k₂, l₂]

        c₁₁₁ = c₁₁₁₁ + (w - w₁) * (c₁₁₁₂ - c₁₁₁₁) / crd[4].Δx
        c₁₂₁ = c₁₂₁₁ + (w - w₁) * (c₁₂₁₂ - c₁₂₁₁) / crd[4].Δx
        c₂₁₁ = c₂₁₁₁ + (w - w₁) * (c₂₁₁₂ - c₂₁₁₁) / crd[4].Δx
        c₂₂₁ = c₂₂₁₁ + (w - w₁) * (c₂₂₁₂ - c₂₂₁₁) / crd[4].Δx
        c₁₁₂ = c₁₁₂₁ + (w - w₁) * (c₁₁₂₂ - c₁₁₂₁) / crd[4].Δx
        c₁₂₂ = c₁₂₂₁ + (w - w₁) * (c₁₂₂₂ - c₁₂₂₁) / crd[4].Δx
        c₂₁₂ = c₂₁₂₁ + (w - w₁) * (c₂₁₂₂ - c₂₁₂₁) / crd[4].Δx
        c₂₂₂ = c₂₂₂₁ + (w - w₁) * (c₂₂₂₂ - c₂₂₂₁) / crd[4].Δx

        c₁₁ = c₁₁₁ + (z - z₁) * (c₁₁₂ - c₁₁₁) / crd[3].Δx
        c₂₁ = c₂₁₁ + (z - z₁) * (c₂₁₂ - c₂₁₁) / crd[3].Δx
        c₁₂ = c₁₂₁ + (z - z₁) * (c₁₂₂ - c₁₂₁) / crd[3].Δx
        c₂₂ = c₂₂₁ + (z - z₁) * (c₂₂₂ - c₂₂₁) / crd[3].Δx

        c₁ = c₁₁ + (y - y₁) * (c₁₂ - c₁₁) / crd[2].Δx
        c₂ = c₂₁ + (y - y₁) * (c₂₂ - c₂₁) / crd[2].Δx

        return c₁ + (x - x₁) * (c₂ - c₁) / crd[1].Δx
    end
end

#####
##### SimpleInterpolation
#####

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

end # module