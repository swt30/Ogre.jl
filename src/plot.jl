# Plotting setup

# This file uses RecipesBase.jl to define plotting behaviour and layout without
# actually loading any plotting backend. To plot a planetary structure, use the
# Plots package and then just call the plot command on a planetary structure. At
# the moment, only planets with atmospheric layers are handled this way - the
# old plotting behaviour is stored in plots-old.jl until it is all transitioned.

using RecipesBase
using LaTeXStrings

# Plotting planetary structures

@recipe function plot(p::PlanetStructure{WithAtmosphere})
    M = vec(mass(p))
    r = vec(radius(p))
    P = vec(pressure(p))
    T = vec(temperature(p))
    τ = vec(opticaldepth(p))

    legend := false
    markershape := :circle
    markersize := 1
    linewidth := 2

    # T-P profile
    # xlims := (300, Inf)
    # ylims := (1e-3, 1e3)
    # yflip := true
    # yscale := :log10
    # (T/K, P/bar)

    # τ-P profile
    xguide := L"Optical depth $τ$"
    yguide := "Pressure P / bar"
    ylims := (1e-5, 1e2)
    xscale := :log10
    yscale := :log10
    (τ, P/bar)

    # # all profiles
    # layout := 4
    # xguide := L"Radius r / R$_⊕$"
    # yguide := hcat(L"Mass m / M$_⊕$",
    #                "Pressure P / bar",
    #                "Temperature T / K",
    #                L"Optical depth $τ$")
    # ylims := [nothing]
    # yscale := [:identity :log10 :identity :log10]
    #
    # x = r/R_earth
    # ys = hcat(M/M_earth, P/1bar, T/K, τ)
    #
    # (x, ys)
end
