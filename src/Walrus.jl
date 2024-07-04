module Walrus

export SimpleInterpolation

export WallStress, WallStressBoundaryConditions

export HomogeneousBodyHeating, PARModelHeating

export Tide, TidalForcing

export WindStress, WindStressBoundaryConditions, LogarithmicNeutralWind

export SurfaceHeatExchange, SurfaceHeatExchangeBoundaryCondition

export WindDrivenStokesDrift, WindDrivenStokesDriftSetup

using Adapt: adapt

import Adapt: adapt_structure

# This is how we will standarise function vs values, I will have to think of a way to deal with discrete vs continuous at some point
struct ReturnValue{FT}
    value :: FT
end

adapt_structure(to, rv::ReturnValue) = ReturnValue(adapt(to, rv.value))

@inline (return_value::ReturnValue)(args...) = return_value.value

display_input(return_value::ReturnValue) = return_value.value
display_input(::Function) = "Function"

include("interpolations.jl")
include("wall_model.jl")
include("radiative_transfer/radiative_transfer.jl")
include("tidal_forcings.jl")
include("wind_stress.jl")
include("surface_heating.jl")
include("wind_driven_stokes.jl")

using .Interpolations
using .WallStressModel
using .RadiativeTransfer
using .TidalForcings
using .WindStressModel
using .SurfaceHeatingModel
using .WindDrivenStokesParameterisation

end # module Walrus
