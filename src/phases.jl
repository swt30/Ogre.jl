module Phases

# Imports
import WaterData
import GeometricalPredicates, VoronoiDelaunay, Dierckx # JLD issues
using BasicUnits
  const bar = BasicUnits.bar
  const cm = BasicUnits.cm
using DataFrames
using Ogre
import GeometricalPredicates: geta, getb, getc, getx, gety
import WaterData: lognorm12, unlognorm, findtriangle, barycoords, Polygon, BoundingBox
using RecipesBase, PlotUtils, LaTeXStrings
using Parameters

# Exports
export findphase, PlanetPhases

# EOS
piecewise = WaterData.load_piecewise_eoses()
fe = piecewise["fe"]
mgsio3 = piecewise["mgsio3"]
h2o = WaterData.load_full_eos()["gridPlusIdeal"]

# Planet structure
# include with default values 70% core, core = 1/3 iron, 2/3 silicate
@with_kw type PlanetModel
  M::Float64 = M_earth
  f_iron::Float64 = 0.7/3
  f_silicate::Float64 = 0.7 - f_iron
  f_water::Float64 = 1 - f_iron - f_silicate
  Tirr::Float64 = 300.
  ɛ::Float64 = 1e-14
  κ = nothing
  γ::Float64 = 0.01
  Psurf::Float64 = 100bar
end

# distinguishing the core & water layers
# here, "core" refers to both iron+silicate
isiron(M, Miron) = (M < Miron)
issilicate(M, Mcore, Miron) = (Miron <= M < Mcore)
iswater(M, Mcore) = (M >= Mcore)
function findphase_cored(M, Mcore, Miron, P, T)
  if iswater(M, Mcore)
    findphase(P, T)
  elseif issilicate(M, Mcore, Miron)
    Silicate
  elseif isiron(M, Mcore)
    Iron
  else
    error("Couldn't get phase")
  end
end
function structure(pm::PlanetModel)
  @unpack M, f_iron, f_silicate, Tirr, ɛ, κ, γ, Psurf = pm
  structure, radius = Ogre.interior(M, f_iron, f_silicate, ɛ, Tirr, κ, γ, Psurf)
  Ts = Ogre.temperature(structure) |> reverse
  Ps = Ogre.pressure(structure) |> reverse
  Rs = Ogre.radius(structure) |> reverse
  Ms = Ogre.mass(structure) |> reverse
  clamp!(Rs, 0, Inf)

  return Ms, Rs, Ps, Ts
end
function radius(pm::PlanetModel)
  Ms, Rs, Ps, Ts = structure(pm)
  return maximum(Rs)
end
function phases(pm::PlanetModel)
  Ms, Rs, Ps, Ts = structure(pm)
  Mcore = pm.M * (pm.f_iron + pm.f_silicate)
  Miron = pm.M * pm.f_iron
  phases = map(Ms, Ps, Ts) do M, P, T
    findphase_cored(M, Mcore, Miron, P, T)
  end
  return PlanetPhases(Rs, phases)
end
function finddensity_cored(M, Mcore, Miron, P, T)
  if iswater(M, Mcore)
    h2o(P, T)
  elseif issilicate(M, Mcore, Miron)
    mgsio3(P, T)
  elseif isiron(M, Mcore)
    fe(P, T)
  else
    error("Couldn't get density")
  end
end
function densities(pm::PlanetModel)
  Ms, Rs, Ps, Ts = structure(pm)
  Mcore = pm.M * (pm.f_iron + pm.f_silicate)
  Miron = pm.M * pm.f_iron
  densities = map(Ms, Ps, Ts) do M, P, T
    finddensity_cored(M, Mcore, Miron, P, T)
  end

  return PlanetInterior(Rs, densities)
end

"Phases of water (or core material)"
@enum Phase Gas Liquid SupercriticalFluid I II III V VI VII VIII X Plasma Superionic Unknown Silicate Iron
"Mappings from phases to plot colors (pretty)"
const phase_colors = Dict(
  Gas => :Azure,
  Liquid => :SkyBlue,
  SupercriticalFluid => :DeepSkyBlue,
  I => :White,
  II => :White,
  III => :White,
  V => :White,
  VI => :MediumPurple,
  VII => :MediumSlateBlue,
  VIII => :DarkSlateBlue,
  X => :MidnightBlue,
  Plasma => :Coral,
  Superionic => :Gold,
  Unknown => :OrangeRed,
  Silicate => "#383232",
  Iron => :Black )

"Mappings from strings to a Phase"
const phasemappings = Dict(
  "idealgas"=>Gas,
  "iapws"=>Unknown, "french"=>Unknown, "fallback"=>Unknown,
  "fluid"=>SupercriticalFluid,
  "I"=>I,
  "II"=>II,
  "III"=>III,
  "V"=>V,
  "VI"=>VI,
  "VII"=>VII, "ice VII"=>VII, "sugimura"=>VII,
  "VIII"=>VIII,
  "X"=>X, "ice X"=>X,
  "superionic"=>Superionic,
  "plasma"=>Plasma )
Base.convert(::Type{Phase}, s::String) = phasemappings[s]

"Critical pressure of water"
Pcrit = WaterData.Pc
"Critical temperature of water"
Tcrit = WaterData.Tc
# regions for water phases
dummy_eos = WaterData.ConstantEOS(1.)
liquid = let
  bounds = WaterData.BoundingBox(0, Inf, 0, Tcrit)
WaterData.BoundedEOS(dummy_eos, bounds)
end
gas = let
bounds = WaterData.BoundingBox(0, Pcrit, Tcrit, Inf)
WaterData.BoundedEOS(dummy_eos, bounds)
end
supercritical = let
bounds = WaterData.BoundingBox(Pcrit, Inf, Tcrit, Inf)
WaterData.BoundedEOS(dummy_eos, bounds)
end

"EOSes with defined regions and associated phases"
immutable EOSCollectionPhases
  eoses::Vector{WaterData.EOS}
  phases::Vector{Phase}
end
eostables = let
  eosdata = WaterData.load_full_eos()
  eosraw = eosdata["raw"]
  eosgrid = eosdata["grid"]
  idealgas = eosdata["gridPlusIdeal"].eoses[1]

  eoses = copy(eosraw.eoses)[[2, 4, 5, 6, 7, 8, 10, 11]]

  eos_strings = ["sugimura", "I", "III", "V", "VI", "VII", "X", "X"]
  eos_phases = [phasemappings[p] for p in eos_strings]

  # add ideal gas
  push!(eoses, idealgas)
  push!(eos_phases, Gas)

  # add other phases
  push!(eoses, liquid)
  push!(eos_phases, Liquid)
  push!(eoses, supercritical)
  push!(eos_phases, SupercriticalFluid)
  push!(eoses, gas)
  push!(eos_phases, Gas)

  EOSCollectionPhases(eoses, eos_phases)
end
function Base.getindex(c::EOSCollectionPhases, P, T)::Phase
  for (eos, phase) in zip(c.eoses, c.phases)
    if (P, T) in eos
      return phase
    end
  end
  return Unknown
end

"An unstructured EOS table with associated phases for each point"
immutable TablePhases <: WaterData.TabularEOS
  P::Vector{Float64}
  T::Vector{Float64}
  ρ::Vector{Float64}
  phase::Vector{String}
  tess::WaterData.Tessellation

  function TablePhases(P, T, ρ, phase)
    u = WaterData.unique_indices(zip(P,T))
    new(P[u], T[u], ρ[u], phase[u], WaterData.get_tessellation(P[u], T[u], uselog=true))
  end
end
function Base.in(P, T, eos::TablePhases)
  if !in(P, T, BoundingBox(eos))
    return false
  else
    xn, yn = lognorm12(P, T, eos)
    return in(xn, yn, eos.tess)
  end
end
function closestpoint(p::TablePhases, P, T)
  Prange = extrema(p.P)
  Trange = extrema(p.T)
  PP = lognorm12(P, Prange...)
  TT = lognorm12(T, Trange...)
  tri = findtriangle(p.tess, PP, TT)
  i = indmax(barycoords(tri, PP, TT))
  closest = if i == 1
    geta(tri)
  elseif i == 2
    getb(tri)
  elseif i == 3
    getc(tri)
  end

  Pclose = unlognorm(getx(closest), Prange...)
  Tclose = unlognorm(gety(closest), Trange...)

  return (Pclose, Tclose)
end
function Base.getindex(p::TablePhases, P, T)::Phase
  if (P, T) in p
    (Pclose, Tclose) = closestpoint(p, P, T)
    for (Pi, Ti, phase) in zip(p.P, p.T, p.phase)
      if Pi ≈ Pclose && Ti ≈ Tclose
        return phase
      end
    end
  end
  return Unknown
end

"The French EOS and its associated phases"
frenchphases = let
  df = readtable(WaterData.config.rawdata * "/French.eos", allowcomments=true)
  P = collect(df[:pressure]) * bar * 1e3
  T = collect(df[:temperature]) * K
  ρ = collect(df[:density]) * g/cm^3
  phase = collect(df[:state])

  TablePhases(P, T, ρ, phase)
end

# run through all phase sources to find the appropriate phase code
function findphase(P, T)::Phase
  thephase = Unknown
  allphases = [frenchphases, eostables]
  for p in allphases
    if p[P, T] != Unknown
      thephase = p[P, T]
      break
    end
  end
  return thephase
end

function change_indices(x)
  where_changed = (x[1:end-1] .!= x[2:end])
  push!(where_changed, true) # ending counts as a change
  find(where_changed)
end

type PlanetPhases
  rs::Vector{Float64}
  ϕs::Vector{Phase}
end

type PlanetInterior
  rs::Vector{Float64}
  xs::Vector{Float64}
end

@recipe function plot(p::PlanetPhases)
  rs, ϕs = p.rs, p.ϕs
  unshift!(ϕs, ϕs[1])
  unshift!(rs, 0)
  colors = Symbol[]
  for (p, c) in phase_colors
    unshift!(rs, 0)
    unshift!(ϕs, p)
    unshift!(colors, c)
  end

  Θs = linspace(0, 2π)
  phasegrid = repeat(ϕs, inner=(1, length(Θs)))
  phasegrid_str = string.(phasegrid)
  Rchanges = rs[change_indices(ϕs)]
  Rp = maximum(rs)
  Rencl = ceil(Rp/R_earth)

  legend := false
  proj := :polar
  border := nothing
  yaxis := nothing
  # colorbar := false
  xaxis := nothing
  ylims --> (0, Rencl)

  # the phase structure
  @series begin
    seriestype := :heatmap
    color := ColorGradient(colors)

    Θs, rs/R_earth, phasegrid_str
  end

  # the borders
  @series begin
    seriestype := :path
    color := :black
    linewidth := 1

    Θs, Rchanges'/R_earth
  end

  # the outer border
  @series begin
    seriestype := :path
    color := :black
    linewidth := 2

    Θs, [Rp/R_earth]
  end

  # the inner grid
  @series begin
    seriestype := :path
    color := :grey
    linestyle := :dot
    linewidth := 0.5

    Θs, collect(1:Rencl-1)'
  end

  # the outer grid
  @series begin
    seriestype := :path
    color := :grey
    linestyle := :dot
    linewidth := 0.5

    Θs, collect(Rencl:10)'
  end
end

end # module Phases
