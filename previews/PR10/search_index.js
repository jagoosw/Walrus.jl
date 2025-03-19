var documenterSearchIndex = {"docs":
[{"location":"appendix/library/#library_api","page":"Library","title":"Library","text":"","category":"section"},{"location":"appendix/library/","page":"Library","title":"Library","text":"Documenting the user interface.","category":"page"},{"location":"appendix/library/#Wall-stress-model","page":"Library","title":"Wall stress model","text":"","category":"section"},{"location":"appendix/library/","page":"Library","title":"Library","text":"Modules = [Walrus.WallStressModel]\nprivate = true","category":"page"},{"location":"appendix/library/#Walrus.WallStressModel","page":"Library","title":"Walrus.WallStressModel","text":"WallStressModel\n\nA wall stress model for LES simulation with default parameters similar to that proposed in Schumann (1975), Hartel (1996), Piomelli et al. (1989), and Taylor and Sarkar (2007).\n\n\n\n\n\n","category":"module"},{"location":"appendix/library/#Walrus.WallStressModel.WallStress-Union{Tuple{}, Tuple{FT}} where FT","page":"Library","title":"Walrus.WallStressModel.WallStress","text":"WallStress(; von_Karman_constant = 0.4,\n             kinematic_viscosity = 1e-6,\n             B = 5.2,\n             precomputed_friction_velocities = false,\n             precompute_speeds = [0:25/100000:25;],\n             grid = nothing,\n             arch = isnothing(grid) ? CPU() : architecture(grid))\n\nReturns a wall stress model for LES simulation with default parameters similar to that proposed in Schumann (1975), Hartel (1996), Piomelli et al. (1989), and Taylor and Sarkar (2007).\n\nFriction velocities will be precomputed at precompute_speeds if  precomputed_friction_velocities is true and grid is provided.\n\nKeyword Arguments\n\nvon_Karman_constant: the von Karman wall stress constant\nkinematic_viscosity: kinematic viscosity of the water above the wall\nB: wall stress constant \nprecomputed_friction_velocities: precompute friction velocities?\nprecompute_speeds: bottom water speeds to precompute friction velocities for,  this should encompas the range of speeds possible in your simulation\ngrid: the grid to precompute the friction velocities for\narch: architecture to adapt precomputed friction velocities for\n\nExample\n\njulia> using Walrus: WallStress\n\njulia> using Oceananigans\n\njulia> wall_stress = WallStress()\n(::WallStress{Float64, Nothing}) (generic function with 1 method)\njulia> boundary_conditions = (u = FieldBoundaryConditions(bottom = FluxBoundaryCondition(wall_stress, discrete_form = true, parameters = Val(:x))),\n                              v = FieldBoundaryConditions(bottom = FluxBoundaryCondition(wall_stress, discrete_form = true, parameters = Val(:y))))\n(u = Oceananigans.FieldBoundaryConditions, with boundary conditions\n├── west: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── east: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── south: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── north: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── bottom: FluxBoundaryCondition: DiscreteBoundaryFunction (::WallStress{Float64, Nothing}) with parameters Val{:x}\n├── top: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n└── immersed: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing), v = Oceananigans.FieldBoundaryConditions, with boundary conditions\n├── west: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── east: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── south: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── north: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── bottom: FluxBoundaryCondition: DiscreteBoundaryFunction (::WallStress{Float64, Nothing}) with parameters Val{:y}\n├── top: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n└── immersed: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing))\n\n\n\n\n\n\n","category":"method"},{"location":"appendix/library/#Walrus.WallStressModel.WallStressBoundaryConditions-Union{Tuple{}, Tuple{FT}} where FT","page":"Library","title":"Walrus.WallStressModel.WallStressBoundaryConditions","text":"WallStressBoundaryConditions(; von_Karman_constant = 0.4,\n                               kinematic_viscosity = 1e-6,\n                               B = 5.2,\n                               precompute_speeds = [0:25/100000:25;],\n                               grid = nothing)\n\nConvenience constructor to setup WallStress boundary conditions.\n\nKeyword Arguments\n\nvon_Karman_constant: the von Karman wall stress constant\nkinematic_viscosity: kinematic viscosity of the water above the wall\nB: wall stress constant \nprecomputed_friction_velocities: precompute friction velocities?\nprecompute_speeds: bottom water speeds to precompute friction velocities for,  this should encompas the range of speeds possible in your simulation\ngrid: the grid to precompute the friction velocities for\n\nExample\n\njulia> using Walrus: WallStressBoundaryConditions\n\njulia> using Oceananigans\n\njulia> stress_boundary_conditions = WallStressBoundaryConditions()\n(u = FluxBoundaryCondition: DiscreteBoundaryFunction (::WallStress{Float64, Nothing}) with parameters Val{:x}, v = FluxBoundaryCondition: DiscreteBoundaryFunction (::WallStress{Float64, Nothing}) with parameters Val{:y})\njulia> boundary_conditions = (u = FieldBoundaryConditions(bottom = stress_boundary_conditions.u),\n                              v = FieldBoundaryConditions(bottom = stress_boundary_conditions.v))\n(u = Oceananigans.FieldBoundaryConditions, with boundary conditions\n├── west: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── east: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── south: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── north: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── bottom: FluxBoundaryCondition: DiscreteBoundaryFunction (::WallStress{Float64, Nothing}) with parameters Val{:x}\n├── top: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n└── immersed: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing), v = Oceananigans.FieldBoundaryConditions, with boundary conditions\n├── west: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── east: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── south: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── north: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── bottom: FluxBoundaryCondition: DiscreteBoundaryFunction (::WallStress{Float64, Nothing}) with parameters Val{:y}\n├── top: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n└── immersed: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing))\n\n\n\n\n\n","category":"method"},{"location":"appendix/library/#Radiative-transfer-models","page":"Library","title":"Radiative transfer models","text":"","category":"section"},{"location":"appendix/library/","page":"Library","title":"Library","text":"Modules = [Walrus.RadiativeTransfer]\nprivate = true","category":"page"},{"location":"appendix/library/#Walrus.RadiativeTransfer","page":"Library","title":"Walrus.RadiativeTransfer","text":"RadiativeTransfer\n\nIncludes models for ratiative transfer through water which can induce body heating\n\n\n\n\n\n","category":"module"},{"location":"appendix/library/#Walrus.RadiativeTransfer.HomogeneousBodyHeating","page":"Library","title":"Walrus.RadiativeTransfer.HomogeneousBodyHeating","text":"HomogeneousBodyHeating\n\nA model for single band light attenuation which heats the water in the form:\n\nI(x y z) = I_0(x y) * expleft(-alpha zright)\n\nwhere I is the radiation intensity and alpha is the attenuation coefficient. This heats the water like fracpartial T(x y z)partial t = fracI(x y z)Ac^prho where A is the area of the cell, c^p is the specific heat capacity and rho is the water density. \n\n\n\n\n\n","category":"type"},{"location":"appendix/library/#Walrus.RadiativeTransfer.HomogeneousBodyHeating-Tuple{}","page":"Library","title":"Walrus.RadiativeTransfer.HomogeneousBodyHeating","text":"HomogeneousBodyHeating(; surface_flux,\n                         water_attenuation_coefficient = 1.8,\n                         water_heat_capacity = 3991.0, # J K⁻¹ kg⁻¹\n                         water_density = 1026.0) # kg m⁻³\n\nCreates a model in which a surface_flux (W / m²) is attenuated by and heats the water. This interacts with Oceananigans as a body forcing.\n\nKeyword Arguments\n\nsurface_flux (required): a function returning the surface radiaiton flux in the form surface_flux(x, y, t) or single value\nwater_attenuation_coefficient: the radiation attenuation coefficient of the water\nwater_heat_capacity: the specific heat capacity of the water\nwater_density: density of the water\n\nExample\n\njulia> using Walrus: HomogeneousBodyHeating\n\njulia> using Oceananigans\n\njulia> grid = RectilinearGrid(size = (128, 128, 128), extent = (1000, 1000, 1000))\n128×128×128 RectilinearGrid{Float64, Oceananigans.Grids.Periodic, Oceananigans.Grids.Periodic, Oceananigans.Grids.Bounded} on Oceananigans.Architectures.CPU with 3×3×3 halo\n├── Periodic x ∈ [0.0, 1000.0)  regularly spaced with Δx=7.8125\n├── Periodic y ∈ [0.0, 1000.0)  regularly spaced with Δy=7.8125\n└── Bounded  z ∈ [-1000.0, 0.0] regularly spaced with Δz=7.8125\njulia> body_heating = HomogeneousBodyHeating(; surface_flux = (x, y, t) -> 100)\n(::HomogeneousBodyHeating{Float64, var\"#1#2\"}) (generic function with 1 method)\n\njulia> model = NonhydrostaticModel(; grid, forcing = (; T = Forcing(body_heating, discrete_form=true)), tracers = :T)\nNonhydrostaticModel{CPU, RectilinearGrid}(time = 0 seconds, iteration = 0)\n├── grid: 128×128×128 RectilinearGrid{Float64, Oceananigans.Grids.Periodic, Oceananigans.Grids.Periodic, Oceananigans.Grids.Bounded} on Oceananigans.Architectures.CPU with 3×3×3 halo\n├── timestepper: QuasiAdamsBashforth2TimeStepper\n├── advection scheme: Centered reconstruction order 2\n├── tracers: T\n├── closure: Nothing\n├── buoyancy: Nothing\n└── coriolis: Nothing\n\n\n\n\n\n","category":"method"},{"location":"appendix/library/#Tidal-forcing-models","page":"Library","title":"Tidal forcing models","text":"","category":"section"},{"location":"appendix/library/","page":"Library","title":"Library","text":"Modules = [Walrus.TidalForcings]\nprivate = true","category":"page"},{"location":"appendix/library/#Walrus.TidalForcings","page":"Library","title":"Walrus.TidalForcings","text":"TidalForcing\n\nProvides quick setup of tidal forcing\n\n\n\n\n\n","category":"module"},{"location":"appendix/library/#Walrus.TidalForcings.Tide","page":"Library","title":"Walrus.TidalForcings.Tide","text":"Tide(; x_amplitude, \n       y_amplitude, \n       period = 12.3782216453hours,\n       nodal_time = 0., \n       x_lag = 0., \n       y_lag = 0.,\n       coriolis = nothing)\n\nSets up a model of tidal forcing with default parameters of an M2 tide.\n\nKeyword Arguments\n\nx_amplitude: the tidal amplitude in the x direction\ny_amplitude: the tidal amplitude in the x direction\nperiod: the tidal period (defaults to that of an M2 tide)\nnodal_time: the time at which peak flow occurs\nx_lag: the phase lag for the tidal component in the x direction\ny_lag: the phase lag for the tidal component in the y direction\ncoriolis: a model for the coriolis parameter \n\nExample\n\njulia> using Walrus: Tide\n\njulia> using Oceananigans\n\njulia> tide = Tide(x_amplitude = 0.1, y_amplitude = 0.)\n(::Tide{Float64, Nothing}) (generic function with 2 methods)\njulia> forcing = (u = Forcing(tide, parameters = Val(:x), discrete_form = true),\n                  v = Forcing(tide, parameters = Val(:y), discrete_form = true))\n(u = DiscreteForcing{Val{:x}}\n├── func: (::Tide{Float64, Nothing}) (generic function with 2 methods)\n└── parameters: Val{:x}(), v = DiscreteForcing{Val{:y}}\n├── func: (::Tide{Float64, Nothing}) (generic function with 2 methods)\n└── parameters: Val{:y}())\n\n\n\n\n\n","category":"type"},{"location":"appendix/library/#Walrus.TidalForcings.TidalForcing-Tuple{}","page":"Library","title":"Walrus.TidalForcings.TidalForcing","text":"TidalForcing(; x_amplitude,\n               y_amplitude,\n               period = 12.3782216453hours,\n               nodal_time = 0.,\n               x_lag = 0.,\n               y_lag = 0.,\n               coriolis = nothing)\n\nA convenience constructor for Tide which returns the forcings pre wrapped.\n\nKeyword Arguments\n\nx_amplitude: the tidal amplitude in the x direction\ny_amplitude: the tidal amplitude in the x direction\nperiod: the tidal period (defaults to that of an M2 tide)\nnodal_time: the time at which peak flow occurs\nx_lag: the phase lag for the tidal component in the x direction\ny_lag: the phase lag for the tidal component in the y direction\ncoriolis: a model for the coriolis parameter \n\nExample\n\njulia> using Walrus: TidalForcing\n\njulia> tidal_forcing = TidalForcing(x_amplitude = 0.1, y_amplitude = 0.)\n(u = DiscreteForcing{Val{:x}}\n├── func: (::Tide{Float64, Nothing}) (generic function with 2 methods)\n└── parameters: Val{:x}(), v = DiscreteForcing{Val{:y}}\n├── func: (::Tide{Float64, Nothing}) (generic function with 2 methods)\n└── parameters: Val{:y}())\n\n\n\n\n\n","category":"method"},{"location":"appendix/library/#Wind-stress-model","page":"Library","title":"Wind stress model","text":"","category":"section"},{"location":"appendix/library/","page":"Library","title":"Library","text":"Modules = [Walrus.WindStressModel]\nprivate = true","category":"page"},{"location":"appendix/library/#Walrus.WindStressModel.LogarithmicNeutralWind-Union{Tuple{}, Tuple{FT}} where FT","page":"Library","title":"Walrus.WindStressModel.LogarithmicNeutralWind","text":"LogarithmicNeutralWind(; monin_obukhov_stability_length = 0.4\n                         charnock_coefficient = 0.014\n                         air_kinematic_viscosity = 1.488e-5\n                         gravity_wave_coefficient = 0.11\n                         gravity = g_Earth,\n\n                         precomputed_roughness_length = false,\n                         precompute_wind_speeds = [0:25/100000:25;],\n                         arch = CPU())\n\nReturns a LogarithmicNeutralWind parameterisation for the surface drag coefficient\n\nC_d is parameterised as,\n\nC_d = left(frackappalogfrac10z_0right)^2\n\nwhere kappa is the Monin‐Obukhov stability length and z_0 is the velocity  roughness length. This is the roughness length scale which logarithmically brings  the relative velocity to zero at the surface, i.e.\n\nU=fracustarkappalogfraczz_0\n\nwhere ustar is the friction velocity. Additionally z_0 is given as,\n\nz_0=bfracnuustar + fraca_cgustar^2\n\nwhere nu is the kinematic viscosity of air and g is the acceleration of gravity.\n\nThis model itterativly solves these equations to find z_0. Alternativly, if the flag  precomputed_roughness_length is set to they are pre computed at precompute_wind_speeds  between which z_0 is then interpolated during run time. Precomputed velocities are  converted to appropriate types for arch (i.e. CPU() or GPU())\n\nThis parameterisaion is described in Smith (1988)\n\n\n\n\n\n","category":"method"},{"location":"appendix/library/#Walrus.WindStressModel.WindStress-Tuple{}","page":"Library","title":"Walrus.WindStressModel.WindStress","text":"WindStress(; reference_wind_speed, \n             reference_wind_direction,\n             drag_coefficient = LogarithmicNeutralWind(), \n             air_density = 1.225, \n             water_density = 1026.)\n\nReturns a wind stress model where the stress is given by,\n\nfractaurho_o = rho_aC_dSU_xy\n\nwhere rho_o is the water density, rho_a is the air density, C_d is the drag coefficient, U_xy are the x and y components of  relative wind speed, and S=sqrtU_x^2+U_y^2.\n\nC_d is calculated from a parameterisation, by default this is a \"log neutral\" wind parameterisation with velocity roughness length parameterisaion like Smith (1988).\n\nIn the default configuration this is the same as described in Fairall et al. (2011).\n\nKeyword Arguments\n\nreference_wind_speed (required): a function returning the (10m neutral) wind speed in the form reference_wind_speed(x, y, t) or single value\nreference_wind_direction (required): a function returning the (10m neutral) wind direction in the form reference_wind_direction(x, y, t) or single value\ndrag_coefficient: the drag coefficient parameterisation\nair_density: air density in kg/m³ \nwater_density: water density in kg/m³ \n\nExample\n\njulia> using Walrus: WindStress\n\njulia> using Oceananigans\n\njulia> reference_wind_speed = 0.1\n0.1\n\njulia> reference_wind_direction = 0.\n0.0\n\njulia> wind_stress = WindStress(; reference_wind_speed, reference_wind_direction)\n(::WindStress{Float64, Float64, LogarithmicNeutralWind{Float64, Nothing}, Float64}) (generic function with 2 methods)\n\njulia> boundary_conditions = (u = FieldBoundaryConditions(top = FluxBoundaryCondition(wind_stress, parameters = Val(:x))),\n                              v = FieldBoundaryConditions(top = FluxBoundaryCondition(wind_stress, parameters = Val(:y))))\n(u = Oceananigans.FieldBoundaryConditions, with boundary conditions\n├── west: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── east: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── south: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── north: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── bottom: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── top: FluxBoundaryCondition: ContinuousBoundaryFunction (::WindStress{Float64, Float64, LogarithmicNeutralWind{Float64, Nothing}, Float64}) at (Nothing, Nothing, Nothing)\n└── immersed: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing), v = Oceananigans.FieldBoundaryConditions, with boundary conditions\n├── west: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── east: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── south: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── north: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── bottom: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── top: FluxBoundaryCondition: ContinuousBoundaryFunction (::WindStress{Float64, Float64, LogarithmicNeutralWind{Float64, Nothing}, Float64}) at (Nothing, Nothing, Nothing)\n└── immersed: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing))\n\n\n\n\n\n\n","category":"method"},{"location":"appendix/library/#Walrus.WindStressModel.WindStressBoundaryConditions-Tuple{}","page":"Library","title":"Walrus.WindStressModel.WindStressBoundaryConditions","text":"WindStressBoundaryConditions(; reference_wind_speed, \n                               reference_wind_direction,\n                               drag_coefficient = LogarithmicNeutralWind(), \n                               air_density = 1.225, \n                               water_density = 1026.)\n\nConvenience constructor to setup WindStress boundary conditions.\n\nKeyword Arguments\n\nreference_wind_speed (required): a function returning the (10m neutral) wind speed in the form reference_wind_speed(x, y, t) or single value\nreference_wind_direction (required): a function returning the (10m neutral) wind direction in the form reference_wind_direction(x, y, t) or single value\ndrag_coefficient: the drag coefficient parameterisation\nair_density: air density in kg/m³ \nwater_density: water density in kg/m³ \n\nExample\n\njulia> using Walrus: WindStressBoundaryConditions\n\njulia> using Oceananigans\n\njulia> wind_stress_boundary_conditions = WindStressBoundaryConditions(; reference_wind_speed = 0.1, reference_wind_direction = 90.)\n(u = FluxBoundaryCondition: DiscreteBoundaryFunction (::WindStress{Float64, Float64, LogarithmicNeutralWind{Float64, Nothing}, Float64}) with parameters Val{:x}, v = FluxBoundaryCondition: DiscreteBoundaryFunction (::WindStress{Float64, Float64, LogarithmicNeutralWind{Float64, Nothing}, Float64}) with parameters Val{:y})\n\njulia> boundary_conditions = (u = FieldBoundaryConditions(top = wind_stress_boundary_conditions.u),\n                              v = FieldBoundaryConditions(top = wind_stress_boundary_conditions.v))\n(u = Oceananigans.FieldBoundaryConditions, with boundary conditions\n├── west: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── east: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── south: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── north: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── bottom: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── top: FluxBoundaryCondition: DiscreteBoundaryFunction (::WindStress{Float64, Float64, LogarithmicNeutralWind{Float64, Nothing}, Float64}) with parameters Val{:x}\n└── immersed: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing), v = Oceananigans.FieldBoundaryConditions, with boundary conditions\n├── west: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── east: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── south: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── north: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── bottom: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── top: FluxBoundaryCondition: DiscreteBoundaryFunction (::WindStress{Float64, Float64, LogarithmicNeutralWind{Float64, Nothing}, Float64}) with parameters Val{:y}\n└── immersed: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing))\n\n\n\n\n\n\n","category":"method"},{"location":"appendix/library/#Walrus.WindStressModel.find_velocity_roughness_length-Tuple{LogarithmicNeutralWind{<:Any, Nothing}, Any, Any, Any}","page":"Library","title":"Walrus.WindStressModel.find_velocity_roughness_length","text":"find_velocity_roughness_length(wind_speed, reference_height, params)\n\nA function that finds the velocity roughness length for the LogarithmicNeutralWind drag coefficient model.\n\nThis will sometimes fail as the function is not well behaved at either low reference heights (it has been tuned for 10m wind), or high (⪆ 20 m/s).\n\n\n\n\n\n","category":"method"},{"location":"appendix/library/#Surface-heat-exchange-model","page":"Library","title":"Surface heat exchange model","text":"","category":"section"},{"location":"appendix/library/","page":"Library","title":"Library","text":"Modules = [Walrus.SurfaceHeatingModel]\nprivate = true","category":"page"},{"location":"appendix/library/#Walrus.SurfaceHeatingModel.SurfaceHeatExchange-Tuple{}","page":"Library","title":"Walrus.SurfaceHeatingModel.SurfaceHeatExchange","text":"SurfaceHeatExchange(; wind_stress,\n                      air_temperature = 18, # °C\n                      latent_heat_vaporisation = EmpiricalLatentHeatVaporisation(),\n                      vapour_pressure = AugustRocheMagnusVapourPressure(),\n                      water_specific_heat_capacity = 3991., # J / K / kg\n                      water_density = 1026., # kg / m³\n                      air_specific_heat_capacity = 1003.5, # J / K / kg\n                      air_density = 1.204, # kg\n                      air_water_mixing_ratio = 0.001, # kg / kg\n                      stephan_boltzman_constant = 5.670374419e-8) # W / K⁴\n\nSpecifies surface heat exhange in the form:\n\nQ = Qᵢᵣ + Qₛ + Qₗ\n\nwhere Qᵢᵣ is the heat flux due to long wave (infra red) radiation, Qₛ is the sensible heat flux, and Qₗ is the latent heat flux. Notably the short wave radiation flux is neglegted here as it is assumed to penatrate far enough into the  water that it is treated separately by HomogeneousBodyHeating (we therefore are also assuming that the short wave penitration is suitably short that it is neglected).\n\nThe short wave term is given by the Stephan-Boltzman equation so:\n\nQᵢᵣ = σ(T⁴ - Tₐ⁴)\n\nwhere σ is the Stephan-Boltzman constant, T is the ocean temperature, and Tₐ is the 2m air temperature. \n\nThe sensible and latent heat flux's are given by the bluk parameterisations described in  Fairall et al. (2011) and are given by:\n\nQₛ = ρₐcₚₐCₕS(T - Tₐ)\n\nand Qₗ = ρₐLₑCₕS(q(T) - qₐ),  where ρₐ is the density of the air, cₚₐ is the specific heat capacity of air,  Cₕ is the heat transfer coefficient, S is the wind speed, Lₑ is the latent heat of vaporizaion parameterised by latent_heat_vaporisation, q(T) is the  saturation vapour pressure of water and is parameterised by vapour_pressure.\n\nThe heat transfer coefficients are given by the same parameterisation as the WindStress. For example for the LogarithmicNeutralWind the transfer coefficient is given as: C_h = frackappalogfrac2z_0frackappalogfrac2z_ot, where z_ot is the sclar roughness parameter given by: z_ot = minleft(115cdot10^-4 55cdot10^-5R_r^-06right), where R_r is the roughness reynolds number given as R_r = fracustar z_0nu where nu is the kinematic viscosity of air.\n\nThe heat flux is then given by:\n\nF = fracQrho_oc_po\n\nwhere rho_o and c_po are the density and specific heat capacity of water.\n\n(Note: we will retain the Oceananigans convention that negative heat flux   at a top boundary increases temeprature).\n\nKeyword Arguments\n\nwind_stress: wind stress model\nair_temperature: the air temperature in °C as a function with signature (x, y, t) or a constant\nlatent_heat_vaporisation: latent heat of vaporisation in J / kg\nvapour_pressure: parameterisation for saturation vapour pressure in water\nwater_specific_heat_capacity: the specific heat capacity of water in J / K / kg\nwater_density: water density in kg / m³\nair_specific_heat_capacit: the specific heat capacity of air in J / K / kg\nair_density: air density in kg / m³\nair_water_mixing_ratio: water content of air in kg / kg\nstephan_boltzman_constant: the Stephan-Boltzman constant in W / K⁴\n\nExample\n\njulia> using Walrus\n\njulia> using Oceananigans\n\njulia> wind_stress = WindStress(; reference_wind_speed = 0., reference_wind_direction = 90.)\n(::WindStress{Float64, Float64, LogarithmicNeutralWind{Float64, Nothing}, Float64}) (generic function with 2 methods)\n\njulia> surface_heat_exchange = SurfaceHeatExchange(; wind_stress)\n(::SurfaceHeatExchange{WindStress{Float64, Float64, LogarithmicNeutralWind{Float64, Nothing}, Float64}, Int64, Walrus.SurfaceHeatingModel.EmpiricalLatentHeatVaporisation{Float64}, Walrus.SurfaceHeatingModel.AugustRocheMagnusVapourPressure{Float64}, Float64}) (generic function with 1 method)\n\njulia> boundary_conditions = (; T = FieldBoundaryConditions(top = FluxBoundaryCondition(surface_heat_exchange, field_dependencies = (:T, :u, :v))))\n(T = Oceananigans.FieldBoundaryConditions, with boundary conditions\n├── west: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── east: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── south: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── north: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── bottom: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing)\n├── top: FluxBoundaryCondition: ContinuousBoundaryFunction (::SurfaceHeatExchange{WindStress{Float64, Float64, LogarithmicNeutralWind{Float64, Nothing}, Float64}, Int64, Walrus.SurfaceHeatingModel.EmpiricalLatentHeatVaporisation{Float64}, Walrus.SurfaceHeatingModel.AugustRocheMagnusVapourPressure{Float64}, Float64}) at (Nothing, Nothing, Nothing)\n└── immersed: DefaultBoundaryCondition (FluxBoundaryCondition: Nothing),)\n\n\n\n\n\n\n","category":"method"},{"location":"appendix/library/#Walrus.SurfaceHeatingModel.SurfaceHeatExchangeBoundaryCondition-Tuple{}","page":"Library","title":"Walrus.SurfaceHeatingModel.SurfaceHeatExchangeBoundaryCondition","text":"SurfaceHeatExchange(; wind_stress,\n                      air_temperature = 18, # °C\n                      latent_heat_vaporisation = EmpiricalLatentHeatVaporisation(),\n                      vapour_pressure = AugustRocheMagnusVapourPressure(),\n                      water_specific_heat_capacity = 3991., # J / K / kg\n                      water_density = 1026., # kg / m³\n                      air_specific_heat_capacity = 1003.5, # J / K / kg\n                      air_density = 1.204, # kg\n                      air_water_mixing_ratio = 0.001, # kg / kg\n                      stephan_boltzman_constant = 5.670374419e-8) # W / K⁴\n\nA convenience constructor returning SurfaceHeatExchange as a boundary condition\n\nKeyword Arguments\n\nwind_stress: wind stress model\nair_temperature: the air temperature in °C as a function with signature (x, y, t) or a constant\nlatent_heat_vaporisation: latent heat of vaporisation in J / kg\nvapour_pressure: parameterisation for saturation vapour pressure in water\nwater_specific_heat_capacity: the specific heat capacity of water in J / K / kg\nwater_density: water density in kg / m³\nair_specific_heat_capacit: the specific heat capacity of air in J / K / kg\nair_density: air density in kg / m³\nair_water_mixing_ratio: water content of air in kg / kg\nstephan_boltzman_constant: the Stephan-Boltzman constant in W / K⁴\n\nExample\n\njulia> using Walrus\n\njulia> wind_stress = WindStress(; reference_wind_speed = 0., reference_wind_direction = 90.)\n(::WindStress{Float64, Float64, LogarithmicNeutralWind{Float64, Nothing}, Float64}) (generic function with 2 methods)\n\njulia> surface_heat_exchange = SurfaceHeatExchangeBoundaryCondition(; wind_stress)\nFluxBoundaryCondition: DiscreteBoundaryFunction with (::SurfaceHeatExchange{WindStress{Float64, Float64, LogarithmicNeutralWind{Float64, Nothing}, Float64}, Int64, Walrus.SurfaceHeatingModel.EmpiricalLatentHeatVaporisation{Float64}, Walrus.SurfaceHeatingModel.AugustRocheMagnusVapourPressure{Float64}, Float64})\n\n\n\n\n\n\n","category":"method"},{"location":"references/#References","page":"References","title":"References","text":"","category":"section"},{"location":"references/","page":"References","title":"References","text":"Fairall, C. W.; Yang, M.; Bariteau, L.; Edson, J. B.; Helmig, D.; McGillis, W.; Pezoa, S.; Hare, J. E.; Huebert, B. and Blomquist, B. (2011). Implementation of the Coupled Ocean-Atmosphere Response Experiment flux algorithm with CO2, dimethyl sulfide, and O3. Journal of Geophysical Research: Oceans 116.\n\n\n\nHartel, C. (1996). Chapter 5 - Turbulent flows: Direct numerical simulation and large-eddy simulation. In: Handbook of Computational Fluid Mechanics, edited by Peyret, R. (Academic Press, London); pp. 283–338.\n\n\n\nPiomelli, U.; Ferziger, J.; Moin, P. and Kim, J. (1989). New approximate boundary conditions for large eddy simulations of wall‐bounded flows. Physics of Fluids A: Fluid Dynamics 1, 1061–1068, arXiv:https://pubs.aip.org/aip/pof/article-pdf/1/6/1061/12366739/1061_1_online.pdf.\n\n\n\nSchumann, U. (1975). Subgrid scale model for finite difference simulations of turbulent flows in plane channels and annuli. Journal of Computational Physics 18, 376–404.\n\n\n\nSmith, S. D. (1988). Coefficients for sea surface wind stress, heat flux, and wind profiles as a function of wind speed and temperature. Journal of Geophysical Research: Oceans 93, 15467–15472.\n\n\n\nTaylor, J. R. and Sarkar, S. (2007). Internal gravity waves generated by a turbulent bottom Ekman layer. Journal of Fluid Mechanics 590, 331–354.\n\n\n\n","category":"page"},{"location":"coming-soon/#Coming-soon","page":"Coming soon","title":"Coming soon","text":"","category":"section"},{"location":"coming-soon/","page":"Coming soon","title":"Coming soon","text":"...","category":"page"},{"location":"appendix/function_index/#Index","page":"Function index","title":"Index","text":"","category":"section"},{"location":"appendix/function_index/","page":"Function index","title":"Function index","text":"","category":"page"},{"location":"#Walrus","page":"Home","title":"Walrus","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Walrus.jl is a Julia package providing various closures models (closure -> seal -> seal 🦭 -> walrus) for ocean flavoured fluid dynamics with Oceananigans.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Currently, it includes:","category":"page"},{"location":"","page":"Home","title":"Home","text":"WallStress: a wall stress model for LES\nTidalForcing: simple body forcing to replicate tides\nRadiativeHeating: a (very) simple radiative transfer model that imparts a body heating\nWindStress: a bulk formulation for wind stress\nSurfaceHeatExchange: a bulk formulation for heat exchange at the ocean surface","category":"page"},{"location":"","page":"Home","title":"Home","text":"Details of the parametrisations used for each model and implementation notes can be found in their docstrings and in the function library.","category":"page"},{"location":"","page":"Home","title":"Home","text":"If you have any questions or suggestions please get in touch through an issue.","category":"page"}]
}
