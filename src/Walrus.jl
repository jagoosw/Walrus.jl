module Walrus

export WallStress, WallStressBoundaryConditions

export HomogeneousBodyHeating

export Tide, TidalForcing

export WindStress, WindStressBoundaryConditions, LogarithmicNeutralWind

export SurfaceHeatExchange, SurfaceHeatExchangeBoundaryCondition

# This is how we will standarise function vs values, I will have to think of a way to deal with discrete vs continuous at some point
struct ReturnValue{FT}
    value :: FT
end

@inline (return_value::ReturnValue)(args...) = return_value.value

display_input(return_value::ReturnValue) = return_value.value
display_input(::Function) = "Function"

include("wall_model.jl")
include("radiative_transfer/radiative_transfer.jl")
include("tidal_forcings.jl")
include("wind_stress.jl")
include("surface_heating.jl")

using .WallStressModel
using .RadiativeTransfer
using .TidalForcings
using .WindStressModel
using .SurfaceHeatingModel

end # module Walrus
