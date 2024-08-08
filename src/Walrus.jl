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

include("get_value.jl")
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
