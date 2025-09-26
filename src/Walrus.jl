module Walrus

export SimpleInterpolation

export WallStress, WallStressBoundaryConditions

export HomogeneousBodyHeating, PARModelHeating

export Tide, Tides

#export WindStress, WindStressBoundaryConditions

#export SurfaceHeatExchange, SurfaceHeatExchangeBoundaryCondition

export WindDrivenStokesDrift, WindDrivenStokesDriftSetup

using Adapt: adapt

import Adapt: adapt_structure

include("get_value.jl")
include("interpolations.jl")
include("wall_model.jl")
include("radiative_transfer/radiative_transfer.jl")
include("tidal_forcings.jl")
include("SurfaceFluxModel/SurfaceFluxModel.jl")
#include("wind_stress.jl")
#include("wind_driven_stokes.jl")
#include("surface_heating.jl")

using .Interpolations
using .WallStressModel
using .RadiativeTransfer
using .TidalForcing
#using .WindStressModel
#using .SurfaceHeatingModel
#using .WindDrivenStokesParameterisation

end # module Walrus
