using Walrus, Oceananigans, OceanBioME, Test

arch = CPU()

include("test_interpolation.jl")
include("test_wind_model.jl")
include("test_wall_model.jl")
include("test_heating.jl")