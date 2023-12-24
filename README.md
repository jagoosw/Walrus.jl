# Walrus 🦭

[![Documentation](https://img.shields.io/badge/documentation-stable%20release-blue?style=flat-square)](https://walrus.jagosw.com/stable/)
[![Documentation](https://img.shields.io/badge/documentation-dev%20release-orange?style=flat-square)](https://walrus.jagosw.com/dev/)
[![Test status badge](https://github.com/jagoosw/Walrus.jl/actions/workflows/test.yml/badge.svg)](https://github.com/jagoosw/Walrus.jl/actions/workflows/test.yml)

``Walrus.jl`` is a Julia package providing various closures models (closure -> seal -> seal 🦭 -> walrus) for ocean flavoured fluid dynamics with [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl/).

Currently, it includes:
- `WallStress`: a wall stress model for LES
- `TidalForcing`: simple body forcing to replicate tides
- `RadiativeHeating`: a (very) simple radiative transfer model that imparts a body heating
- `WindStress`: a bulk formulation for wind stress
- `SurfaceHeatExchange`: a bulk formulation for heat exchange at the ocean surface

Details of the parametrisations used for each model and implementation notes can be found in their docstrings and in [the documentation](https://walrus.jagosw.com/stable/appendix/library/).

If you have any questions or suggestions please get in touch through an [issue](https://github.com/jagoosw/Walrus.jl/issues).

If you use this package in your work please contact me so that I can expedite the process of making it citable.