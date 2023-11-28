module Walrus

export WallStress, WallStressBoundaryConditions

export BodyHeating

include("wall_model.jl")
include("radiative_transfer/radiative_transfer.jl")
include("tidal_forcing.jl")

using .WallStressModel
using .RadiativeTransfer

end # module Walrus
