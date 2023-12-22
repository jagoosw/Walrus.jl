using Walrus, Oceananigans, Documenter, Test

include("test_windstress.jl")
include("test_wall_model.jl")

doctest(Walrus)