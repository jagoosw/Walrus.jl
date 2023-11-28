module Walrus

export WallStress, WallStressBoundaryConditions

include("wall_model.jl")
include("tidal_forcing.jl")

using .WallStressModel

end # module Walrus
