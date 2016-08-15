# Plotting setup

# This file uses PlotRecipes.jl to define plotting behaviour and layout without
# actually loading any plotting backend. To plot a planetary structure, use the
# Plots package and then just call the plot command on a planetary structure. At
# the moment, only planets with atmospheric layers are handled this way - the
# old plotting behaviour is stored in plots-old.jl until it is all transitioned.

using PlotRecipes
using LaTeXStrings

# Bring in the EOS so we can plot density as a variable
import WaterData
h2o_full = WaterData.load_full_eos()["grid"]
h2o_ideal = WaterData.load_functional_eoses()["misc"]["ideal_gas"]
function ρ_patched(P, T)
    if P < P_rad_max
        return h2o_full(P, T)
    else
        return h2o_ideal(P, T)
    end
end

# Plotting planetary structures
@recipe function plot(p::PlanetStructure{WithAtmosphere})
    M = vec(mass(p))
    r = vec(radius(p))
    P = vec(pressure(p))
    T = vec(temperature(p))
    τ = vec(opticaldepth(p))
    ρ = map(ρ_patched, P, T)

    legend := false
    markershape := :xcross
    markersize := 3
    linewidth := 2

    # T-P profile
    # xlims := (300, Inf)
    # ylims := (1e-3, 1e3)
    # yflip := true
    # yscale := :log10
    # (T/K, P/bar)

    # τ-P profile
    # xguide := L"Optical depth $τ$"
    # yguide := "Pressure P / bar"
    # ylims := (1e-3, 1e3)
    # (τ, P/bar)

    # all profiles
    layout := 5
    xguide := L"Radius r / R$_⊕$"
    yguide := hcat(L"Mass m / M$_⊕$",
                   "Pressure P / bar",
                   "Temperature T / K",
                   L"Optical depth $τ$",
                   L"Density $ρ$ / kg/m$^3$")
    ylims := [nothing]
    yscale := [:identity :log10 :identity :log10 :log10]

    x = r/R_earth
    ys = hcat(M/M_earth, P/1bar, T/K, τ, ρ/(kg/m^3))

    (x, ys)
end
