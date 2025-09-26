using Walrus: normalise_surface_function, get_value

struct PrescribedAtmosphericState{WS, WD, AT, MR, DL, FT}
                  wind_speed :: WS 
              wind_direction :: WD
                 temperature :: AT
      air_water_mixing_ratio :: MR
        downwelling_longwave :: DL

                 wind_height :: FT
          temperature_height :: FT
end

function PrescribedAtmosphericState(; wind_speed, 
                                      wind_direction,
                                      temperature,
                                      air_water_mixing_ratio = 0.006, # kg / kgs
                                      downwelling_longwave = EmpiricalDownwellingLongwave(),

                                      wind_height = 10,
                                      temperature_height = 10)
        
    wind_speed = normalise_surface_function(wind_speed)
    wind_direction = normalise_surface_function(wind_direction)
    temperature = normalise_surface_function(temperature)
    air_water_mixing_ratio = normalise_surface_function(air_water_mixing_ratio)

    return PrescribedAtmosphericState(wind_speed, wind_direction,
                                      temperature, air_water_mixing_ratio,
                                      downwelling_longwave,
                                      wind_height, temperature_height)
end

@inline function x_wind(atmosphere::PrescribedAtmosphericState, i, j, grid, clock, model_fields)
    U = get_value(atmosphere.wind_speed, i, j, grid, clock)
    θ = get_value(atmosphere.wind_direction, i, j, grid, clock)

    return - U * sind(θ)
end

@inline relative_x_wind(atmosphere::PrescribedAtmosphericState, i, j, grid, clock, model_fields) =
    (@inbounds x_wind(atmosphere, i, j, grid, clock, model_fields) - model_fields.u[i, j, grid.Nz])

@inline function y_wind(atmosphere::PrescribedAtmosphericState, i, j, grid, clock, model_fields)
    U = get_value(atmosphere.wind_speed, i, j, grid, clock)
    θ = get_value(atmosphere.wind_direction, i, j, grid, clock)

    return - U * cosd(θ)
end

@inline relative_y_wind(atmosphere::PrescribedAtmosphericState, i, j, grid, clock, model_fields) =
   (@inbounds y_wind(atmosphere, i, j, grid, clock, model_fields) -  model_fields.v[i, j, grid.Nz])

@inline wind_speed(atmosphere::PrescribedAtmosphericState, i, j, grid, clock, model_fields) =
    get_value(atmosphere.wind_speed, i, j, grid, clock)

@inline function relative_wind_speed(atmosphere::PrescribedAtmosphericState, i, j, grid, clock, model_fields)
    u = @inbounds model_fields.u[i, j, grid.Nz]
    v = @inbounds model_fields.v[i, j, grid.Nz]

    uʷ = x_wind(atmosphere, i, j, grid, clock, model_fields)
    vʷ = y_wind(atmosphere, i, j, grid, clock, model_fields)

    return √((uʷ - u)^2 + (vʷ - v)^2)
end

@inline temperature(atmosphere::PrescribedAtmosphericState, i, j, grid, clock, model_fields) =
    get_value(atmosphere.temperature, i, j, grid, clock)

@inline air_water_mixing_ratio(atmosphere::PrescribedAtmosphericState, i, j, grid, clock, model_fields) =
    get_value(atmosphere.air_water_mixing_ratio, i, j, grid, clock)

@inline downwelling_longwave(atmosphere::PrescribedAtmosphericState, i, j, grid, clock, model_fields, boundary) =
    atmosphere.downwelling_longwave(i, j, grid, clock, model_fields, atmosphere, boundary)

@inline velocity_reference_height(atmosphere::PrescribedAtmosphericState, i, j, grid, clock, model_fields) = 
    atmosphere.wind_height

@inline temperature_reference_height(atmosphere::PrescribedAtmosphericState, i, j, grid, clock, model_fields) = 
    atmosphere.temperature_height
    
#####
##### empirical downwelling longwave
#####

struct EmpiricalDownwellingLongwave{FT, CF} # Brunt, 1932 / Yang et al., 2023 Atmos. Chem. Phys.
                 a :: FT
                 b :: FT
                 α :: FT
                 β :: FT
                 γ :: FT
                 δ :: FT
                 ζ :: FT
    cloud_fraction :: CF

    function EmpiricalDownwellingLongwave(; a::FT = 0.599,
                                            b::FT = 0.053,
                                            α::FT = 0.178,
                                            β::FT = 0.339,
                                            γ::FT = 0.075,
                                            δ::FT = 0.395,
                                            ζ::FT = 0.253,
                                            cloud_fraction = 0.3) where FT

        cloud_fraction = normalise_surface_function(cloud_fraction)

        CF = typeof(cloud_fraction)

        return new{FT, CF}(a, b, α, β, γ, δ, ζ, cloud_fraction)
    end

  EmpiricalDownwellingLongwave(a::FT, b::FT, α::FT, β::FT, γ::FT, δ::FT, ζ::FT, cloud_fraction::CF) where {FT, CF} =
      new{FT, CF}(a, b, α, β, γ, δ, ζ, cloud_fraction)
end

adapt_structure(to, ed::EmpiricalDownwellingLongwave) = 
  EmpiricalDownwellingLongwave(adapt(to, ed.a),
                               adapt(to, ed.b),
                               adapt(to, ed.α),
                               adapt(to, ed.β),
                               adapt(to, ed.γ),
                               adapt(to, ed.δ),
                               adapt(to, ed.ζ), 
                               adapt(to, ed.cloud_fraction))
  

@inline function (ed::EmpiricalDownwellingLongwave)(i, j, grid, clock, model_fields, atmosphere, boundary)
    a = ed.a
    b = ed.b

    α = ed.α
    β = ed.β
    γ = ed.γ
    δ = ed.δ
    ζ = ed.ζ

    σ = boundary.stephan_boltzman_constant

    T = temperature(atmosphere, i, j, grid, clock, model_fields)

    q  = air_water_mixing_ratio(atmosphere, i, j, grid, clock, model_fields)

    q′ = boundary.vapour_pressure

    e = q *  q′.p / (q′.a + (1 - q′.a) * q)

    N = get_value(ed.cloud_fraction, i, j, grid, clock)

    eₛ = water_vapour_pressure(q′, T)
    RH = 100 * e / eₛ

    return ((a + b * √e) * (1 - α * N^β) + γ * N^δ * RH^ζ) * σ * (T + 273.15)^4
end