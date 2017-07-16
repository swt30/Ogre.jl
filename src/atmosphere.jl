# Atmospheric model which sits on top of the interior structure model

import ForwardDiff: derivative
import Optim: optimize, OnceDifferentiable, Fminbox, minimizer

# New types

"This system includes both temperature and optical depth details"
struct WithAtmosphere <: ModelComplexity end

"An opacity in units of m^2 / kg"
abstract type MassOpacity end

"A power law opacity of the form κ = C P^α T^β (in SI units)"
struct PowerLawOpacity{Op<:Opacity} <: MassOpacity
    C::Op
    α::Float64
    β::Float64
end
ConstantOpacity(C) = PowerLawOpacity(C, 0.0, 0.0)

abstract type OpacityRatio end

struct ConstantOpacityRatio <: OpacityRatio
    γ::Float64
end
(c::ConstantOpacityRatio)(P, T) = c.γ

struct VariableOpacityRatio <: OpacityRatio
    γ::Function
end
(v::VariableOpacityRatio)(P, T) = v.γ(P, T)

"The optical depth gradient w.r.t. mass, dτ/dm = -κ/4πr²"
struct OpticalDepthGradient <: StructureEquation
    κ::PowerLawOpacity
    is_radiative::Function
end

"Temperature gradient for a radiative two-stream atmosphere"
mutable struct TwoStreamTemperatureGradient{T1<:Temperature, T2<:Temperature} <: StructureEquation
    κ::PowerLawOpacity
    Tint::T1
    Tirr::T2
    γ::OpacityRatio
end

"Temperature gradient for a combined atmosphere and interior model"
struct CombinedTemperatureGradient <: StructureEquation
    atmosphere_gradient::TwoStreamTemperatureGradient
    interior_gradient::TemperatureGradient
    is_radiative::Function
end


"A planet with temperature dependence and an atmosphere"
mutable struct AtmospherePlanet <: PlanetSystem{WithAtmosphere}
    M::Float64
    structure_equations::EquationSet
    boundary_values::BoundaryValues{WithAtmosphere}
    solution_grid::Vector{Float64}
    radius_search_bracket::Vector{Float64}
    refine_surface!::Nullable{Function}
end

"Planetary structure that holds temperature and opacity info"
mutable struct AtmospherePlanetStructure <: PlanetStructure{WithAtmosphere}
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
mutable struct AtmosphereValues{Ma<:Mass,Ra<:Radius,
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

function (odg::OpticalDepthGradient)(vs::ValueSet)
    if isphysical(vs) && odg.is_radiative(vs)
        r = radius(vs)
        P = pressure(vs)
        T = temperature(vs)
        κ = odg.κ(P, T)

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
        κ = tst.κ(P, T)
        Tint = tst.Tint
        Tirr = tst.Tirr
        γ = tst.γ(P, T)
        dTdτ = temperature_profile_deriv(τ, Tint, Tirr, γ)
        dTdm = dTdτ * dτdm(κ, r)
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

const water_opacity_r_vis = kurosaki_opacity(2.20, 1.0, -0.4)
const water_opacity_r_th = kurosaki_opacity(3.07e2, 0.9, -4.0)
const water_opacity_p_vis = kurosaki_opacity(1.94e4, 0.01, 1.0)
const water_opacity_p_th = kurosaki_opacity(4.15e5, 0.01, -1.1)

const water_opacity = water_opacity_r_th
const α = water_opacity_r_th.α
const β = water_opacity_r_th.β
const C = water_opacity_r_th.C * m^2 / kg
const water_opacity_ratio = VariableOpacityRatio((P,T) -> water_opacity_r_vis(P, T) / water_opacity_r_th(P, T))

is_radiative_default(vs::ValueSet) = (pressure(vs) < P_rad_max)

function scale_height(T, g)
    H = (R_H2O * T / g)
    return H
end

function surface_opticaldepth(M, R, T, γ)
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

function surface_temperature(τsurf, Tint, Tirr, γ)
    temperature_profile(τsurf, Tint/K, Tirr/K, γ) * K
end

function surface_T_and_τ(Tint, Tirr, M, R, γ)
    τguess = 0.01
    Tguess = Tirr
    τf(T) = surface_opticaldepth(M, R, T, γ)
    Tf(τ) = surface_temperature(τ, Tint, Tirr, γ)
    lsq(T, τ) = (Tf(τ)/K - T)^2 + (τf(T*K) - τ)^2
    lsq(Tτ) = lsq(Tτ...)
    lbounds = [0., 0.]
    ubounds = [Inf, Inf]
    result = optimize(OnceDifferentiable(lsq), [Tguess/K, τguess],
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
                          r_bracket=defaults.R_bracket, κ_const=nothing,
                          γ_const=nothing, Psurf=P_rad_max, refine_surface=nothing)

    # set up structural equations
    masscontinuity = MassContinuity(eos)
    adiabatic_gradient = TemperatureGradient(eos, Cₚ)

    is_variable_opacity = κ_const == nothing
    is_variable_opacity_ratio = γ_const == nothing
    κ = (is_variable_opacity ? water_opacity : ConstantOpacity(κ_const))
    γ = (is_variable_opacity_ratio ? water_opacity_ratio : ConstantOpacityRatio(γ_const))
    changed_Pradmax = Psurf != P_rad_max
    is_radiative = if changed_Pradmax
        vs -> pressure(vs) < Psurf
    else
        is_radiative_default
    end
    wateratm_gradient = TwoStreamTemperatureGradient(κ, Tint_initial, Tirr, γ)
    opticaldepth_gradient = OpticalDepthGradient(κ, is_radiative)
    combined_temperature = CombinedTemperatureGradient(wateratm_gradient,
                                                       adiabatic_gradient,
                                                       is_radiative)

    structure = EquationSet([masscontinuity, pressurebalance,
                             combined_temperature, opticaldepth_gradient])

    AtmospherePlanet(M, structure, bvs, grid, r_bracket, refine_surface)
end

function mass_grid(M, f, Npoints)
    minusM = [0.]
    append!(minusM, logspace(-12, log10(f), Npoints÷2))
    append!(minusM, linspace(f, 1, Npoints÷2))
    grid = M * (1 - minusM)
end

function mass_grid_linear(M, f, Npoints)
    grid = linspace(M, 0, Npoints)
end
