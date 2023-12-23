using Walrus, Oceananigans, Documenter, Test

include("test_wind_model.jl")
include("test_wall_model.jl")
include("test_heating.jl")

doctest(Walrus)