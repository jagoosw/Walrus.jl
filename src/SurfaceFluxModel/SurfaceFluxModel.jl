module SurfaceFluxModel

export OceanAtmosphereBoundaryConditions, PrescribedAtmosphericState, SimilarityTheoryInterface

using Adapt

import Adapt: adapt_structure
import Oceananigans.BoundaryConditions: update_boundary_condition!

include("similarity_theory_coefficients.jl")
include("prescribed_atmospheric_state.jl")
include("boundary.jl")
include("precomputed_coefficients.jl")

end # moduleg