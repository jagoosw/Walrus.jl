using Walrus, Oceananigans, OceanBioME, Test, CUDA

arch = CPU()

@info "Testing on $arch"

include("test_interpolation.jl")
include("test_get_value.jl")
include("test_surface_flux.jl")
include("test_body_heating.jl")
include("test_wall_model.jl")
include("test_heating.jl")
include("test_wind_driven_stokes.jl")