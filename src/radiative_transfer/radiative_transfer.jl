"""
RadiativeTransfer

Includes models for ratiative transfer through water which can induce body heating
"""
module RadiativeTransfer

export BodyHeating

using KernelAbstractions

using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: regularize_field_boundary_conditions, fill_halo_regions!, ValueBoundaryCondition, FieldBoundaryConditions
using Oceananigans.Fields: CenterField
using Oceananigans.Grids: node, znode
using Oceananigans.Utils: launch!

using KernelAbstractions.Extras: @unroll

import Oceananigans.Biogeochemistry: update_biogeochemical_state!, update_tendencies!, AbstractBiogeochemistry

import Base: summary, show

RadiationField(grid, surface_flux) = CenterField(grid; 
                                                 boundary_conditions = 
                                                    regularize_field_boundary_conditions(
                                                        FieldBoundaryConditions(top = ValueBoundaryCondition(surface_flux)),
                                                                                grid, 
                                                                                :radiaton))
include("body_heating.jl")
end # module