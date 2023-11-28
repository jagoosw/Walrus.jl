module Walrus

export WallStress, WallStressBoundaryConditions

export BodyHeating

export Tide, TidalForcing

using Adapt

import Adapt: adapt_structure

include("wall_model.jl")
include("radiative_transfer/radiative_transfer.jl")
include("tidal_forcings.jl")

using .WallStressModel
using .RadiativeTransfer
using .TidalForcings

end # module Walrus
