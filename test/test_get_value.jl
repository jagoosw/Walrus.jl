using Oceananigans.Grids: xnode, ynode

using Walrus: get_value, normalise_surface_function

number_val = 100

continuous_val(x, y, t) = x * y * t
continuous_val(::Nothing, y, t) = y * t

two_d_field_val(x, y) = continuous_val(x, y, 100)
two_d_field_val(y)    = continuous_val(1, y, 100)

three_d_field_val(x, y, z) = continuous_val(x, y, 100)
three_d_field_val(x, y)    = continuous_val(1, y, 100)
three_d_field_val(y)       = continuous_val(1, y, 100)

function discrete_val(i, j, grid, clock)

    x = xnode(i, grid, Center())
    y = ynode(j, grid, Center())

    t = clock.time

    return continuous_val(x, y, t)
end

@testset "Get value" begin
   for grid in [RectilinearGrid(size = (1, 1, 1), extent = (2, 2, 2)),
                RectilinearGrid(topology = (Periodic, Periodic, Flat), size = (1, 1), extent = (2, 2)),
                RectilinearGrid(topology = (Flat, Periodic, Flat), size = (1, ), extent = (2))]

        two_d_field = Field{Center, Center, Nothing}(grid; indices = (:, :, 1))
        three_d_field = Field{Center, Center, Center}(grid)

        set!(two_d_field,   two_d_field_val)
        set!(three_d_field, three_d_field_val)

        clock = Clock(; time = eltype(grid)(100))

        possible_values = [normalise_surface_function.([number_val, continuous_val, two_d_field, three_d_field])..., normalise_surface_function(discrete_val, discrete_form = true)]

        for value in possible_values
            @test get_value(value, 1, 1, grid, clock) == 100
        end
   end
end