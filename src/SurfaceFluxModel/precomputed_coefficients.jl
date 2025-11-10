using KernelAbstractions: @kernel, @index, synchronize
using Oceananigans.Architectures: architecture, on_architecture, device
using Walrus.Interpolations: SimpleInterpolation

struct PrecomputedInterfaceCoefficients{UM, DI, HI, DC, HC}
                         underlying_model :: UM

             drag_coefficient_interpolant :: DI
    heat_exchange_coefficient_interpolant :: HI

                         drag_coefficient :: DC
                heat_exchange_coefficient :: HC
end

function PrecomputedInterfaceCoefficients(grid;
                                          max_iterations = 1000,
                                          underlying_model = SimilarityTheoryInterface(nothing; max_iterations),
                                          wind_speed_range = 0:1:3,#0:1.0:30,
                                          air_temperature_range = 0:1:3,#-40:1.0:40,
                                          mixing_ratio_range = 0:1e-3:3e-3,#0:1e-4:1e-2,
                                          water_temperature_range = 0:1:3,#-5:1.0:40,
                                          velocity_reference_height = 10,
                                          temperature_reference_height = 2,
                                          allow_non_regular_ranges = false)

    arch = architecture(grid)

    FT = eltype(grid)

    ((isregular(wind_speed_range) & 
      isregular(air_temperature_range) & 
      isregular(mixing_ratio_range) & 
      isregular(water_temperature_range)) | allow_non_regular_ranges) ||
        @error "Non-regular ranges are not currently supported for interpolation"

    drag_coefficient_array = on_architecture(arch, zeros(FT, length(wind_speed_range), 
                                                             length(air_temperature_range), 
                                                             length(mixing_ratio_range),
                                                             length(water_temperature_range)))
    heat_exchange_coefficient_array = on_architecture(arch, zeros(FT, length(wind_speed_range), 
                                                                      length(air_temperature_range), 
                                                                      length(mixing_ratio_range),
                                                                      length(water_temperature_range)))

    @info "Precomputing interface coefficients"
    precompute_coefficients! = _precompute_coefficients!(device(arch))

    precompute_coefficients!(underlying_model, 
                             FT.(wind_speed_range), 
                             FT.(air_temperature_range), 
                             FT.(mixing_ratio_range), 
                             FT.(water_temperature_range), 
                             FT.(velocity_reference_height), 
                             FT.(temperature_reference_height),
                             FT.(drag_coefficient_array),
                             heat_exchange_coefficient_array,
                             ndrange = length(drag_coefficient_array))

    synchronize(device(arch))

    drag_coefficient_interpolant = 
        SimpleInterpolation((wind_speed_range, 
                             air_temperature_range,
                             mixing_ratio_range,
                             water_temperature_range),
                             drag_coefficient_array,
                             arch)

    heat_exchange_coefficient_interpolant = 
        SimpleInterpolation((wind_speed_range, 
                             air_temperature_range,
                             mixing_ratio_range,
                             water_temperature_range),
                             heat_exchange_coefficient_array,
                             arch)

    drag_coefficient = Field{Center, Center, Nothing}(grid; indices = (:, :, 1))
    heat_exchange_coefficient = Field{Center, Center, Nothing}(grid; indices = (:, :, 1))

    set!(drag_coefficient, sqrt(1e-3))
    set!(heat_exchange_coefficient, sqrt(1e-3))
    
    return PrecomputedInterfaceCoefficients(underlying_model, 
                                            drag_coefficient_interpolant, 
                                            heat_exchange_coefficient_interpolant, 
                                            drag_coefficient, 
                                            heat_exchange_coefficient)
end

isregular(::Union{<:StepRange, <:StepRangeLen}) = true
isregular(array) = all(diff(array) .== diff(array)[1]) 

@kernel function _precompute_coefficients!(interface::SimilarityTheoryInterface, Ur, θr, wr, Tr, zᵤ, zₜ, Cdr, Chr)
    lin_idx = @index(Global, NTuple)

    i, j, k, l = @inbounds Tuple(CartesianIndices(Cdr)[lin_idx])

    FT = eltype(Ur)

    U = @inbounds Ur[i]
    θ = @inbounds θr[j]
    w = @inbounds wr[k]
    T = @inbounds Tr[l]

    u′, T′ = sqrt(FT(1e-3)), sqrt(FT(1e-3))

    u′₋, T′₋ = FT(Inf), FT(Inf)

    iters = 0
    
    @inbounds while ((abs(u′ - u′₋) > 1e-8) | (abs(T′ - T′₋) > 1e-8)) & (iters <= interface.max_iterations)
        u′₋ = u′
        T′₋ = T′

        next_step = itterate_scaling_values(u′, T′, U, θ, T, w, zᵤ, zₜ, interface)

        u′ = next_step.u′
        T′ = next_step.T′

        iters += 1
    end

    #((abs(u′ - u′₋) > 1e-8) | (abs(T′ - T′₋) > 1e-8)) && (@warn "Did not converge with $u′, $T′, $U, $θ, $T, $w, $(abs(u′ - u′₋)), $(abs(T′ - T′₋))"; u′ = NaN; T′ = NaN)

    Cd = @inbounds u′^2 / (U^2 + eps(0.0))
    Ch = @inbounds - T′ * u′ / (T - θ + eps(0.0)) / (U + eps(0.0))

    @inbounds Cdr[lin_idx] = Cd#min(convert(FT, 1/10), ifelse(U == 0, zero(FT), Cd))
    @inbounds Chr[lin_idx] = Ch#min(convert(FT, 1/10), ifelse(isfinite(Ch), Ch, FT(1e-3)))
end

@kernel function _compute_coefficients!(interface::PrecomputedInterfaceCoefficients, grid, clock, model_fields, atmosphere)
    i, j = @index(Global, NTuple)

    U = relative_wind_speed(atmosphere, i, j, grid, clock, model_fields)
    θ = temperature(atmosphere, i, j, grid, clock, model_fields)
    w = air_water_mixing_ratio(atmosphere, i, j, grid, clock, model_fields)
    T = @inbounds model_fields.T[i, j, grid.Nz]

    interface.drag_coefficient[i, j, 1] = interface.drag_coefficient_interpolant(U, θ, w, T)
    interface.heat_exchange_coefficient[i, j, 1] = interface.heat_exchange_coefficient_interpolant(U, θ, w, T)
end