# Useful physical constants

export G, M_earth, R_earth, σB


"Gravitational constant"
const G = 6.67384e-11 * m^3/(kg*s^2)
"Mass of the Earth, M_⊕"
const M_earth = 5.972e24 * kg
"Radius of the Earth, R_⊕"
const R_earth = 6.3781e6 * m
"""The transition pressure between the radiative and convective layer.
Above this pressure, the convective treatment is used."""
const P_rad_max = 100 * bar
"Constant for an isotropic atmosphere (day-side redistribution)"
const μ_isotropic = 1/√3
"Optical depth of the photosphere"
const τ_photosphere = 1
"Radius of the Sun"
const Rsun = 695700 * km
"Temperature of the Sun"
const Tsun = 5780 * K
"The Stefan-Boltzmann constant"
const σB = 5.670367e-8 * W/(m^2*K^4)
"Molecular weight of water"
const MW_H2O = 18 * g/mol
"Ideal gas constant"
const R_ideal = 8.3144598 * J/(K*mol)
"Specific gas constant for water"
const R_H2O = R_ideal / MW_H2O

# Simple equations using these constants

"Surface gravity"
surface_gravity(M, R) = (G*M/(R^2))
"Surface gravity of Earth"
const g_earth = surface_gravity(M_earth, R_earth)
"Luminosity of a mass uniformly producing energy"
luminosity(M_core, ɛ) = M_core * ɛ
"Flux at a given radius"
flux(L, R) = L / (4π*R^2)
"Effective temperature of a black body"
Teff_blackbody(F) = (F/σB) ^ (1/4)
"Equilibrium temperature for a planet"
function Teq_blackbody(T_star, R_star, semimajor_axis, albedo)
    Ts = T_star
    Rs = R_star
    a = semimajor_axis
    α = albedo

    Ts * √(Rs/2a) * (1 - α)^(1/4)
end
"Surface temperature for a planet based on internal energy"
function Tsurf_from_heat(M_planet, R_planet, core_fraction, heat_per_unit_mass)
    Mp = M_planet
    Rp = R_planet
    Mcore = Mp * core_fraction

    L = luminosity(Mcore, heat_per_unit_mass)
    F = flux(L, Rp)

    Teff = Teff_blackbody(F)
end
