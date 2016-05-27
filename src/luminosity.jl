using Roots: fzero

"Core luminosity from the core mass and energy generation rate"
function luminosity(core_mass, energy_generation_per_mass)
    core_mass * energy_generation_per_mass
end

"Surface flux from an object's luminosity and radius"
function flux(luminosity, radius)
    r = radius
    luminosity / (4π*r^2)
end

"Effective temperature of a blackbody with a given surface flux"
function Teff_blackbody(flux)
    (flux/σB)^(1/4)
end

"Equilibrium temperature for a planet"
function Teq_blackbody(Tstar, Rstar, semimajor_axis, albedo)
    Tstar * sqrt(Rstar/(2*semimajor_axis)) * (1 - albedo)^(1/4)
end

"Total surface temperature for an illuminated and internally heated planet"
Ttotal(Teff, Teq) = (Teff^4 + Teq^4)^(1/4)

"Right-hand side of the photosphere pressure equation, where κ = κ(P, T)"
Prhs(g, κ, P, T) = τ * g / κ(P, T)

"Left-hand side of the photosphere pressure equation"
Plhs(P) = P

"The function to solve to obtain pressure at the photosphere"
Psolve(g, κ, P, T) = Prhs(g, κ, P, T) - Plhs(P)

"Collision-induced Rosseland mean opacity for water at the photosphere"
function κ(P, T)
    # P in Pa
    P_bars = P/10000
    # T in K
    T_1000K = T/1000

    return 2.20 * P_bars^1 * T_100K^(-0.4) * 0.1 # cm^2 / g -> m^2 / kg
end

"Pressure at the photosphere"
function Pphot(g, T)
    fzero(P -> Psolve(g, κ, P, T))
end
