# Atmospheric model which sits on top of the interior structure model

import ForwardDiff: derivative
import Optim: optimize, DifferentiableFunction, Fminbox, minimizer

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
type TwoStreamTemperatureGradient{T1<:Temperature, T2<:Temperature} <: StructureEquation
    κ::PowerLawOpacity
    Tint::T1
    Tirr::T2
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
    return (T⁴ ^ (1/4))::Temperature
end
"Derivative of the angle-dependent atmospheric temperature profile"
function temperature_profile_deriv(τ::OpticalDepth,
                                   T_int::Temperature, T_irr::Temperature,
                                   γ::Dimensionless, μ::Dimensionless=μ_isotropic)
    T = τ -> temperature_profile(τ, T_int, T_irr, γ, μ)
    dTdτ = derivative(T, τ)
end

# Evaluating structure

function (κ::PowerLawOpacity)(P::Pressure, T::Temperature)
    C = κ.C
    α = κ.α
    β = κ.β

    return (C * (P/1Pa)^α * (T/1K)^β)::Opacity
end

"The optical-depth--radius relation for a radiative layer"
function dτdm(κ::Dimensionless, r::Distance)
    return -κ/(4π*r^2)
end

set_isothermal_atm(b::Bool) = (global ISOTHERMAL_ATM = b)
set_fixed_opacity(b::Bool) = (global FIXED_OPACITY = b)
set_κ_const(κ) = (global κ_const = float(κ))
set_P_rad_max(P) = (global P_rad_max = float(P))

const ISOTHERMAL_ATM = false
const FIXED_OPACITY = true
const κ_const = 30.

function (odg::OpticalDepthGradient)(vs::ValueSet)
    if isphysical(vs) && odg.is_radiative(vs)
        r = radius(vs)
        P = pressure(vs)
        T = temperature(vs)
        if FIXED_OPACITY
            κ = κ_const
        else
            κ = odg.κ(P, T)
        end

        return dτdm(κ, r)
    else
        return zero(Float64)
    end
end

function (tst::TwoStreamTemperatureGradient)(vs::ValueSet)
    dTdm = zero(Float64)

    if isphysical(vs)
        r = radius(vs)
        P = pressure(vs)
        T = temperature(vs)
        τ = opticaldepth(vs)
        if FIXED_OPACITY
            κ = κ_const
        else
            κ = tst.κ(P, T)
        end
        Tint = tst.Tint
        Tirr = tst.Tirr
        γ = tst.γ
        if !ISOTHERMAL_ATM
            dTdτ = temperature_profile_deriv(τ, Tint, Tirr, γ)
            dTdm = dTdτ * dτdm(κ, r)
        end
    end

    return dTdm
end

function (ctg::CombinedTemperatureGradient)(vs::ValueSet)
    dTdm = zero(Float64)
    atmosphere = ctg.atmosphere_gradient
    interior = ctg.interior_gradient
    is_radiative = ctg.is_radiative

    dTdm = is_radiative(vs) ? atmosphere(vs) : interior(vs)

    return dTdm
end

# FIXME: this is a workaround for julia issue #14919
# (the macro approach ran into segfault issues)
for structEq in (MassContinuity, PressureBalance, TemperatureGradient,
                 CombinedTemperatureGradient, TwoStreamTemperatureGradient,
                 OpticalDepthGradient)
    @eval begin
        let MRP = Ogre.MassRadiusPressure, PV = Ogre.PhysicalValues, AV=Ogre.AtmosphereValues
            (eq::$structEq)(m, y::Vector) = eq(m, y...)
            (eq::$structEq)(m, r, P) = eq(MRP(m, r, P))
            (eq::$structEq)(m, r, P, T) = eq(PV(m, r, P, T))
            (eq::$structEq)(m, r, P, T, τ) = eq(AV(m, r, P, T, τ))
        end
    end
end

function (tg::TemperatureGradient)(av::AtmosphereValues)
    m = mass(av)
    r = radius(av)
    P = pressure(av)
    T = temperature(av)

    tg(m, r, P, T)
end
(eos::MassPiecewiseEOS)(vs::AtmosphereValues) = extracteos(eos, vs)(vs)

# Other features of the atmospheric model

"""Make a power law opacity from the functions in Kurosaki et al.
These functions have the form κ = D * (P / 1 bar)^α * (T / 1000 K)^β cm^2/g.
This function returns a power law opacity κ(P, T), accepting P/Pa and T/K."""
function kurosaki_opacity(D, α, β)
    C = D * cm^2/g * (1/1bar)^α * (1/1000K)^β
    PowerLawOpacity(C, α, β)
end

const water_opacity_r = kurosaki_opacity(2.20, 1.0, -0.4)
# const water_opacity_r_vis = kurosaki_opacity(2.20, 1.0, -0.4)
# const water_opacity_r_th = kurosaki_opacity(3.07e2, 0.9, -4.0)
# const water_opacity_p_vis = kurosaki_opacity(1.94e4, 0.01, 1.0)
# const water_opacity_p_th = kurosaki_opacity(4.15e5, 0.01, -1.1)
const α = water_opacity_r.α
const β = water_opacity_r.β
const C = water_opacity_r.C * m^2 / kg
const γ = 1.

is_radiative(vs::ValueSet) = (pressure(vs) < P_rad_max)
const opticaldepth_gradient = OpticalDepthGradient(water_opacity_r, is_radiative)

function scale_height(T, g)
    H = (R_H2O * T / g)
    return H
end

function surface_opticaldepth(M, R, T)
    g = surface_gravity(M, R)
    H = scale_height(T, g)
    τ = (1/γ) * sqrt(H / (2π*(α+1)*R))
    return τ
end

function surface_pressure(M, R, T, τ)
    # this function is written in dimensionless form so we divide out units
    let G = G/(m^3/kg/s^2), M = M/kg, R = R/m, C = C/(m^2/kg), T = T/K
        num = G * M * (α+1) * τ
        den = R^2 * C * T^β
        P = (num/den)^(1/(α+1))
        return P * Pa
    end
end

function surface_temperature(τsurf, Tint, Tirr)
    temperature_profile(τsurf, Tint/K, Tirr/K, γ) * K
end

function surface_T_and_τ(Tint, Tirr, M, R)
    τguess = 0.01
    Tguess = Tirr
    τf(T) = surface_opticaldepth(M, R, T)
    Tf(τ) = surface_temperature(τ, Tint, Tirr)
    lsq(T, τ) = (Tf(τ)/K - T)^2 + (τf(T*K) - τ)^2
    lsq(Tτ) = lsq(Tτ...)
    lbounds = [0., 0.]
    ubounds = [Inf, Inf]
    result = optimize(DifferentiableFunction(lsq), [Tguess/K, τguess],
                      lbounds, ubounds, Fminbox())
    T = minimizer(result)[1] * K
    τ = minimizer(result)[2]
    return (T, τ)
end

function T_at_P(ps::PlanetStructure{WithAtmosphere}, P)
    Ps = pressure(ps)
    i = findfirst(x -> x >= P, Ps)
    return temperature(ps)[i]
end

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

# FIXME: this is a workaround for julia issue #14919
macro addEOSCall(eos)
    esc(
        quote
            let P = Ogre.pressure, T = Ogre.temperature,
                MRP = Ogre.MassRadiusPressure, PV = Ogre.PhysicalValues,
                AV = Ogre.AtmosphereValues
                (e::$eos)(vs::MRP) = e(P(vs))
                (e::$eos)(vs::PV) = e(P(vs), T(vs))
                (e::$eos)(vs::AV) = e(P(vs), T(vs))
            end
        end
    )
end

for eos in (TFD, BME3, BME4, Vinet, PolytropicEOS, WaterData.OutOfDomainEOS, BoundedEOS, IAPWS, MGDPressureEOS, PressurePiecewiseEOS, LineEOS, GridEOS, StitchedEOS)
    @addEOSCall eos
end

"Do-everything constructor for a watery planet with a watery atmosphere"
function AtmospherePlanet(M::Mass, eos::EOS, Cₚ::HeatCapacity,
                          bvs::BoundaryValues{WithAtmosphere},
                          Tint_initial::Temperature, Tirr::Temperature,
                          grid=linspace(M, 0, defaults.total_points),
                          r_bracket=defaults.R_bracket, refine_surface=nothing)

    wateratm_gradient = TwoStreamTemperatureGradient(water_opacity_r,
                                                     Tint_initial, Tirr, γ)
    masscontinuity = MassContinuity(eos)
    adiabatic_gradient = TemperatureGradient(eos, Cₚ)
    combined_temperature = CombinedTemperatureGradient(wateratm_gradient,
                                                       adiabatic_gradient,
                                                       is_radiative)

    structure = EquationSet([masscontinuity, pressurebalance,
                             combined_temperature, opticaldepth_gradient])

    AtmospherePlanet(M, structure, bvs, grid, r_bracket, refine_surface)
end

module Heating

export interior
import Ogre, Dierckx, WaterData
import Ogre: M_earth, R_earth
using BasicUnits

# equations of state
simple_eoses = WaterData.load_piecewise_eoses()
const fe = simple_eoses["fe"]
const mgsio3 = simple_eoses["mgsio3"]
const h2o = WaterData.load_full_eos()["grid"]
const Cₚ = WaterData.load_heat_capacity()["heatcap_h2o"]

# integrator settings
const Rbracket = [0, 10] * R_earth
const Rguess = mean(Rbracket)
const Npoints = 500

# make a planet with full heating treatment
function interior(M, f, ɛ, Tirr)
    γ = Ogre.γ

    # Guess initial photospheric parameters
    Tint_guess = Ogre.Tsurf_from_heat(M_earth, R_earth, f, ɛ)
    Tτ_guess = Ogre.surface_T_and_τ(Tint_guess, Tirr, M_earth, R_earth)
    Tphot_guess, τphot_guess = Tτ_guess
    Pphot_guess = Ogre.surface_pressure(M_earth, R_earth,
                                        Tphot_guess, τphot_guess)

    # Define core structure
    massfracs = [1/3, 2/3]
    eoses = WaterData.EOS[fe, mgsio3]
    if f < 1
        massfracs = massfracs * f
        push!(massfracs, 1-f)
        push!(eoses, h2o)
    end
    @assert sum(massfracs) ≈ 1
    planet_eos = Ogre.MassPiecewiseEOS(eoses, M/kg, massfracs)

    # surface boundary conditions
    log10f = log10(f)
    minusM = [0.]
    append!(minusM, logspace(-13, log10f, (2*Npoints÷3)))
    append!(minusM, linspace(f, 1, Npoints÷3))
    grid = M * (1 - minusM)
    bvs = Ogre.ValueSet(M/kg, Rguess/m, Pphot_guess/Pa,
                        Tphot_guess/K, τphot_guess)

    # make the surface refinement function
    refine_surface! = (sys) -> bvs_update!(sys, f, ɛ)

    # make and solve the planet structure
    planet = Ogre.AtmospherePlanet(
        M, planet_eos, Cₚ, bvs,
        Tint_guess, Tirr, grid, Rbracket, refine_surface!)
    Ogre.find_structure_and_radius!(planet)
end

function bvs_update!(sys, f, ɛ)
    bvs = sys.boundary_values
    atmosphere_eqn = Ogre.atm_gradient(sys)

    # get current values
    M = bvs.m
    R = bvs.r
    P = bvs.P
    T = bvs.T
    τ = bvs.τ
    Tirr = atmosphere_eqn.Tirr
    Tint = atmosphere_eqn.Tint

    # update internal temperature
    newTint = Ogre.Tsurf_from_heat(M, R, f, ɛ)
    atmosphere_eqn.Tint = newTint

    # update surface temperature and optical depth
    newT, newτ = Ogre.surface_T_and_τ(newTint, Tirr, M, R)
    bvs.T = newT
    bvs.τ = newτ

    # update surface pressure
    newP = Ogre.surface_pressure(M, R, newT, newτ)
    bvs.P = newP
end

end # module Heating

const interior = Heating.interior
