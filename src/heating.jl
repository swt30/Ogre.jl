module Heating

export interior
import Ogre, Dierckx, WaterData
import Ogre: M_earth, R_earth, defaults, water_opacity_ratio
using BasicUnits

# equations of state
simple_eoses = WaterData.load_piecewise_eoses()
const fe = simple_eoses["fe"]
const mgsio3 = simple_eoses["mgsio3"]
const h2o = WaterData.load_full_eos()["gridPlusIdeal"]
const Cₚ = WaterData.load_heat_capacity()["heatcap_h2o"]

import WaterData.istempdependent
WaterData.istempdependent(::WaterData.PressurePiecewiseEOS) = true

# integrator settings
const Rbracket = defaults.R_bracket
const Rguess = mean(Rbracket)
const Npoints = defaults.total_points

# make a planet with full heating treatment
function interior(M, firon, fsilicate, ɛ, Tirr, κ, γ, Psurf)
    fc = firon + fsilicate
    γ_guess = (γ == nothing) ? water_opacity_ratio(1bar, 300K) : γ
    # Guess initial photospheric parameters
    Tint_guess = Ogre.Tsurf_from_heat(M_earth, R_earth, fc, ɛ)
    Tτ_guess = Ogre.surface_T_and_τ(Tint_guess, Tirr, M_earth, R_earth, γ_guess)
    Tphot_guess, τphot_guess = Tτ_guess
    Pphot_guess = Ogre.surface_pressure(M_earth, R_earth,
                                        Tphot_guess, τphot_guess)

    # Define core structure
    fc_iron = (firon > 0) ? firon/fc : 0
    fc_silicate = (fsilicate > 0) ? fsilicate/fc : 0
    massfracs = [fc_iron, fc_silicate]
    eoses = WaterData.EOS[fe, mgsio3]
    if fc < 1
        massfracs = massfracs * fc
        push!(massfracs, 1-fc)
        push!(eoses, h2o)
    end
    @assert sum(massfracs) ≈ 1
    planet_eos = Ogre.MassPiecewiseEOS(eoses, M/kg, massfracs)

    # surface boundary conditions
    grid = Ogre.mass_grid(M, fc, Npoints)
    bvs = Ogre.ValueSet(M/kg, Rguess/m, Pphot_guess/Pa,
                        Tphot_guess/K, τphot_guess)

    # make the surface refinement function
    refine_surface! = (sys -> bvs_update!(sys, fc, ɛ))

    # make and solve the planet structure
    planet = Ogre.AtmospherePlanet(
        M, planet_eos, Cₚ, bvs,
        Tint_guess, Tirr, grid, Rbracket,
        κ, γ, Psurf, refine_surface!)
    Ogre.find_structure_and_radius!(planet)
end

function bvs_update!(sys, fc, ɛ)
    bvs = sys.boundary_values
    atmosphere_eqn = Ogre.atm_gradient(sys)

    # get current values
    M = bvs.m
    R = bvs.r
    P = bvs.P
    T = bvs.T
    τ = bvs.τ
    γ_surf = atmosphere_eqn.γ(P, T)
    Tirr = atmosphere_eqn.Tirr
    Tint = atmosphere_eqn.Tint

    # update internal temperature
    newTint = Ogre.Tsurf_from_heat(M, R, fc, ɛ)
    atmosphere_eqn.Tint = newTint

    # update surface temperature and optical depth
    newT, newτ = Ogre.surface_T_and_τ(newTint, Tirr, M, R, γ_surf)
    bvs.T = newT
    bvs.τ = newτ

    # update surface pressure
    newP = Ogre.surface_pressure(M, R, newT, newτ)
    bvs.P = newP
end

end # module Heating

const interior = Heating.interior
