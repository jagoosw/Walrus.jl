"""
RadiativeTransfer

Includes models for ratiative transfer through water which can induce body heating
"""
module RadiativeTransfer

export HomogeneousBodyHeating

using KernelAbstractions

using Adapt: adapt
using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: regularize_field_boundary_conditions, fill_halo_regions!, ValueBoundaryCondition, FieldBoundaryConditions
using Oceananigans.Fields: CenterField
using Oceananigans.Forcings: Forcing
using Oceananigans.Grids: node, znode, Center, Face
using Oceananigans.Utils: launch!
using Oceananigans.Operators: Azᶠᶠᶜ

using KernelAbstractions.Extras: @unroll

import Oceananigans.Biogeochemistry: update_biogeochemical_state!, update_tendencies!, AbstractBiogeochemistry

import Adapt: adapt_structure
import Base: summary, show

RadiationField(grid, surface_flux) = CenterField(grid; 
                                                 boundary_conditions = 
                                                    regularize_field_boundary_conditions(
                                                        FieldBoundaryConditions(top = ValueBoundaryCondition(surface_flux)),
                                                                                grid, 
                                                                                :radiaton))
include("homogeneous_body_heating.jl")
end # module