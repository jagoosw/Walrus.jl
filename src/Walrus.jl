module Walrus

export WallStress, WallStressBoundaryConditions

export HomogeneousBodyHeating

export Tide, TidalForcing

using Adapt

include("wall_model.jl")
include("radiative_transfer/radiative_transfer.jl")
include("tidal_forcings.jl")

# This is how we will standarise function vs values, I will have to think of a way to deal with discrete vs continuous at some point
struct ReturnValue{FT}
    value :: FT
end

@inline (return_value::ReturnValue)(args...) = return_value.value

using .WallStressModel
using .RadiativeTransfer
using .TidalForcings

end # module Walrus
