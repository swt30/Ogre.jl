# Atmospheric model which sits on top of the interior structure model

import ForwardDiff: derivative

# Constants

"Constant for an isotropic atmosphere (day-side redistribution)"
const μ_isotropic = 1/√3
"Optical depth of the photosphere"
const τ_photosphere = 1

# New types

"This system includes both temperature and optical depth details"
immutable WithAtmosphere <: ModelComplexity; end

"An opacity in units of m^2 / kg"
abstract MassOpacity

"A power law opacity of the form κ = C P^α T^β (in SI units)"
immutable PowerLawOpacity{Op<:Opacity} <: MassOpacity
    C::Op
    α::Float64
    β::Float64
end

"The optical depth gradient w.r.t. mass, dτ/dm = -κ/4πr²"
immutable OpticalDepthGradient <: StructureEquation
    κ::PowerLawOpacity
    is_radiative::Function
end

"Temperature gradient for a radiative two-stream atmosphere"
type TwoStreamTemperatureGradient{Te<:Temperature} <: StructureEquation
    κ::PowerLawOpacity
    Tint::Te
    Tirr::Te
    γ::Float64
end

"Temperature gradient for a combined atmosphere and interior model"
immutable CombinedTemperatureGradient <: StructureEquation
    atmosphere_gradient::TwoStreamTemperatureGradient
    interior_gradient::TemperatureGradient
    is_radiative::Function
end

"A planet with temperature dependence and an atmosphere"
type AtmospherePlanet <: PlanetSystem{WithAtmosphere}
    M::Float64
    structure_equations::EquationSet
    boundary_values::BoundaryValues{WithAtmosphere}
    solution_grid::Vector{Float64}
    radius_search_bracket::Vector{Float64}
    refine_surface!::Nullable{Function}
end

"Planetary structure that holds temperature and opacity info"
type AtmospherePlanetStructure <: PlanetStructure{WithAtmosphere}
    data::Matrix{Float64}

    function AtmospherePlanetStructure(data::Matrix{Float64})
        if size(data)[1] == 5
            new(data)
        else
            error("Expected a size (5,n) array: got $(size(data))")
        end
    end
end
PlanetStructure(m, r, P, T, τ) = AtmospherePlanetStructure(hcat(m, r, P, T, τ)')

"Holds physical values of mass, radius, pressure, temperature, opacity"
type AtmosphereValues{Ma<:Mass,Ra<:Radius,
                      Pr<:Pressure,Te<:Temperature,
                      Op<:Opacity} <: ValueSet{WithAtmosphere}
    m::Ma
    r::Ra
    P::Pr
    T::Te
    τ::Op
end
ValueSet(m, r, P, T, τ) = AtmosphereValues(m, r, P, T, τ)

# Temperature profiles

"Angle-dependent atmospheric temperature profile"
function temperature_profile(τ::OpticalDepth,
                             T_int::Temperature, T_irr::Temperature,
                             γ::Dimensionless, μ::Dimensionless=μ_isotropic)
    internal_term = 3/4 * T_int^4 * (2/3 + τ)
    irradiation_term = 3/4 * T_irr^4 * μ * (2/3 + μ/γ + (γ/3μ - μ/γ)*exp(-γ*τ/μ))

    T⁴ = internal_term + irradiation_term
    return (T⁴ ^ (1//4))::Temperature
end
"Derivative of the angle-dependent atmospheric temperature profile"
function temperature_profile_deriv(τ::OpticalDepth,
                                   T_int::Temperature, T_irr::Temperature,
                                   γ::Dimensionless, μ::Dimensionless=μ_isotropic)
    T(τ) = temperature_profile(τ, T_int, T_irr, γ, μ)
    dTdτ = derivative(T, τ)
end

# Evaluating structure

function Base.call(κ::PowerLawOpacity, P::Pressure, T::Temperature)
    C = κ.C
    α = κ.α
    β = κ.β

    return (C * (P/1Pa)^α * (T/1K)^β)::Opacity
end

"The optical-depth--radius relation for a radiative layer"
function dτdm(κ::Dimensionless, r::Distance)
    -κ/(4π*r^2)
end

"A constant opacity for testing purposes"
const κ_const = 0.1

function Base.call(odg::OpticalDepthGradient, vs::ValueSet)
    if isphysical(vs) && odg.is_radiative(vs)
        r = radius(vs)
        P = pressure(vs)
        T = temperature(vs)
        # κ = odg.κ(P, T)
        κ = κ_const

        return dτdm(κ, r)
    else
        return zero(Float64)
    end
end

function Base.call(tst::TwoStreamTemperatureGradient, vs::ValueSet)
    if isphysical(vs)
        r = radius(vs)
        P = pressure(vs)
        T = temperature(vs)
        τ = opticaldepth(vs)
        # κ = tst.κ(P, T)
        κ = κ_const
        Tint = tst.Tint
        Tirr = tst.Tirr
        γ = tst.γ

        dTdτ = temperature_profile_deriv(τ, Tint, Tirr, γ)
        dTdm = dTdτ * dτdm(κ, r)
        return dTdm
    else
        return zero(Float64)
    end
end

function Base.call(ctg::CombinedTemperatureGradient, vs::ValueSet)
    atmosphere = ctg.atmosphere_gradient
    interior = ctg.interior_gradient
    is_radiative = ctg.is_radiative

    if is_radiative(vs)
        atmosphere(vs)
    else
        interior(vs)
    end
end

Base.call(e::StructureEquation, m, r, P, T, τ) = e(AtmosphereValues(m, r, P, T, τ))
function Base.call(tg::TemperatureGradient, av::AtmosphereValues)
    m = mass(av)
    r = radius(av)
    P = pressure(av)
    T = temperature(av)

    tg(m, r, P, T)
end
Base.call(eos::MassPiecewiseEOS, vs::AtmosphereValues) = extracteos(eos, vs)(vs)

# Other additions to Ogre methods

function blank_structure(sys::PlanetSystem{WithAtmosphere})
    n = npoints(sys)
    m = nvars(WithAtmosphere)
    AtmospherePlanetStructure(fill(NaN, m, n))
end

atm_gradient(sys::AtmospherePlanet) = sys.structure_equations[3].atmosphere_gradient

function refine_surface!(system::AtmospherePlanet)
    if !isnull(system.refine_surface!)
        do_refine! = get(system.refine_surface!)
        do_refine!(system)
    end
    # if !isnull(system.refine_atmosphere_grid!)
    #     do_refine_grid! = get(system.refine_atmosphere_grid!)
    #     refine_range = 1:10
    #     do_refine_grid!(system.solution_grid, refine_range)
    # end
end

nvars(::Type{WithAtmosphere}) = 5 # M, R, P, T, τ
opticaldepth(av::AtmosphereValues) = av.τ
temperature(av::AtmosphereValues) = av.T
opticaldepth(ps::PlanetStructure{WithAtmosphere}) = ps.data[5, :]
temperature(ps::PlanetStructure{WithAtmosphere}) = ps.data[4, :]
nonmass(ps::PlanetStructure{WithAtmosphere}) = ps.data[2:5, :]
nonmass(av::AtmosphereValues) = [radius(av), pressure(av), temperature(av),
                                 opticaldepth(av)]
setnonmass!(ps::PlanetStructure{WithAtmosphere}, i, rPTτ) = (ps.data[2:5, i] = rPTτ; nothing)
function density(ps::PlanetStructure{WithAtmosphere}, T, sys::PlanetSystem)
    eos = sys.structure_equations[1].eos
    map(eos, pressure(ps), temperature(ps))
end
centre(ps::PlanetStructure{WithAtmosphere}) = AtmosphereValues(ps.data[:, end]...)
surface(ps::PlanetStructure{WithAtmosphere}) = AtmosphereValues(ps.data[:, 1]...)

function Base.show(io::IO, atm::AtmosphereValues)
    m, r, P, T, τ = atm.m, atm.r, atm.P, atm.T, atm.τ
    show(io, "$(m/M_earth) M⊕, $(r/R_earth) R⊕, $(P/Pa) Pa, $(T/K) K, τ=$τ")
end

function isphysical(vs::AtmosphereValues)
    (mass(vs) > 0 && radius(vs) > 0 && pressure(vs) > 0 &&
     temperature(vs) > 0 && opticaldepth(vs) > 0)
end
let specific_eoses = (TFD, BME3, BME4, Vinet, PolytropicEOS, WaterData.OutOfDomainEOS),
    modifier_eoses = (BoundedEOS, WaterData.InverseFunctionalEOS, PressurePiecewiseEOS),
    general_eoses = (LineEOS, GridEOS, EOS, )

    for T in chain(specific_eoses, modifier_eoses, general_eoses)
        Base.call(eos::T, vs::AtmosphereValues) = eos(pressure(vs), temperature(vs))
    end
end

"""
                            Make a power law opacity from the functions in Kurosaki et al.
                            These functions have the form κ = D * (P / 1 bar)^α * (T / 1000 K)^β cm^2/g.
                            This function returns a power law opacity κ(P, T), accepting P/Pa and T/K.
                            """
function kurosaki_opacity(D, α, β)
    C = D * cm^2/g * (1/1bar)^α * (1/1000K)^β
    PowerLawOpacity(C, α, β)
end

# The water opacities
const water_opacity_r_vis = kurosaki_opacity(2.20, 1.0, -0.4)
const water_opacity_r_th = kurosaki_opacity(3.07e2, 0.9, -4.0)
const water_opacity_p_vis = kurosaki_opacity(1.94e4, 0.01, 1.0)
const water_opacity_p_th = kurosaki_opacity(4.15e5, 0.01, -1.1)

# Default values
const Tint = 400. * K
const Tirr = 300. * K
const γ = 1.

is_radiative(vs::ValueSet) = (pressure(vs) < P_rad_max)
const opticaldepth_gradient = OpticalDepthGradient(water_opacity_r_th, is_radiative)
const wateratm_gradient = TwoStreamTemperatureGradient(water_opacity_r_th, Tint, Tirr, γ)

"Do-everything constructor for a watery planet with a watery atmosphere"
function AtmospherePlanet(M::Mass, eos::EOS, Cₚ::HeatCapacity,
                          bvs::BoundaryValues{WithAtmosphere},
                          grid=linspace(M, 0, defaults.total_points),
                          r_bracket=defaults.R_bracket, refine_surface=nothing)

    masscontinuity = MassContinuity(eos)
    adiabatic_gradient = TemperatureGradient(eos, Cₚ)
    combined_temperature = CombinedTemperatureGradient(wateratm_gradient,
                                                       adiabatic_gradient,
                                                       is_radiative)

    structure = EquationSet([masscontinuity, pressurebalance,
                             combined_temperature, opticaldepth_gradient])

    AtmospherePlanet(M, structure, bvs, grid, r_bracket, refine_surface)
end
