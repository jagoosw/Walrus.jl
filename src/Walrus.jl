module Walrus

export WallStress, WallStressBoundaryConditions

export HomogeneousBodyHeating

export Tide, TidalForcing

using Adapt

include("wall_model.jl")
include("radiative_transfer/radiative_transfer.jl")
include("tidal_forcings.jl")

using .WallStressModel
using .RadiativeTransfer
using .TidalForcings

end # module Walrus
