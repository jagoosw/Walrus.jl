using Oceananigans.Fields: AbstractField, Center
using Oceananigans.Grids: xnode, ynode

# fallback
@inline get_value(f, args...) = f

struct ContinuousSurfaceFunction{F}
    func :: F
end

@inline function get_value(f::ContinuousSurfaceFunction, i, j, grid, clock, args...)
    t = clock.time

    x = xnode(i, j, grid.Nz, grid, Center(), Center(), Center())
    y = ynode(i, j, grid.Nz, grid, Center(), Center(), Center())

    return f.func(x, y, t, args...)
end

@inline get_value(f::ContinuousSurfaceFunction, ::Nothing, ::Nothing, t, args...) = f.func(0, 0, t)

struct DiscreteSurfaceFuncton{F}
    func :: F
end

@inline get_value(f::DiscreteSurfaceFuncton, i, j, grid, clock, args...) = f.func(i, j, grid, clock, args...)

# GPU compatible, mainly for fields
@inline get_value(f::AbstractArray{<:Any, 2}, i, j, grid, clock, args...) = @inbounds f[i, j]
@inline get_value(f::AbstractArray{<:Any, 3}, i, j, grid, clock, args...) = @inbounds f[i, j, grid.Nz]

# fallback
normalise_surface_function(f; kwargs...) = f

normalise_surface_function(f::Function; discrete_form = false) = discrete_form ? DiscreteSurfaceFuncton(f) : ContinuousSurfaceFunction(f)
