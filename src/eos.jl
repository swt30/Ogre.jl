# EOS.JL
# Equation of state handling

# Required packages
using Dierckx, Roots

# EOS types
#------------------------------------------------------------------------------

"Equation of State type"
abstract EOS{mc<:ModelComplexity} <: Equation
"Equation of state consisting of just one function"
abstract SingleEOS{mc<:ModelComplexity} <: EOS{mc}
"Equation of state that is subdivided into multiple pieces"
abstract PiecewiseEOS{mc<:ModelComplexity} <: EOS{mc}
"Equation of state that can be calculated directly."
abstract SimpleEOS{mc<:ModelComplexity} <: SingleEOS{mc}

""" An EOS depending solely on pressure.
    
    * `equation`: Function ρ=f(P)
    * `fullname`: Name of the EOS (for printing and plots) """
immutable PressureEOS <: SimpleEOS{NoTemp}
    equation::Function
    fullname::String
end
SimpleEOS(::Type{NoTemp}, f::Function, name::String) = PressureEOS(f, name)


""" An EOS depending on both pressure and temperature

    * `equation`: Function ρ=f(P, T)
    * `fullname`: Name of the EOS (for printing and plots) """
immutable PressureTempEOS <: SimpleEOS{WithTemp}
    equation::Function
    fullname::String
end
SimpleEOS(::Type{WithTemp}, f::Function, name::String) = PressureTempEOS(f, name)


""" Equation of state which inverts some function over a density range

    * `equation`: Function P=f(ρ)
    * `a`, `b`: Density range to invert over
    * `fullname`: Name of the EOS (for printing and plots) """
immutable InvPressureEOS{T<:Real} <: SingleEOS{NoTemp}
    equation::Function
    a::T
    b::T
    fullname::String
end

""" Equation of state which is piecewise in the mass coordinate

    * `equations`: Vector of `SingleEOS` for each piece.
    * `transition_values`: Vector defining the edge of each piece. Evaluating
      past the ends of this vector will just evaluate the nearest EOS.
    * `fullname`: Name of the EOS (for printing and plotting) """
immutable MassPiecewiseEOS{T<:Real, E<:SingleEOS{NoTemp}} <: PiecewiseEOS{NoTemp}
    equations::Vector{E}
    transition_values::Vector{T}
    fullname::String
end
function MassPiecewiseEOS(equations, transition_values)
    fullname = join([n.fullname for n in equations], " & ")
    MassPiecewiseEOS(equations, transition_values, fullname)
end
function MassPiecewiseEOS(equations, M, mass_fractions)
    layer_edges = [0, cumsum(mass_fractions)] .* M
    MassPiecewiseEOS(equations, layer_edges)
end

""" Equation of state which is piecewise in the pressure coordinate

    * `equations`: Vector of `SingleEOS` for each piece.
    * `transition_values`: Vector defining the edge of each piece. Evaluating
      past the ends of this vector will just evaluate the nearest EOS.
    * `fullname`: Name of the EOS (for printing and plotting) """
immutable PressurePiecewiseEOS{T<:Real, E<:SingleEOS} <: PiecewiseEOS
    equations::Vector{E}
    transition_values::Vector{T}
    fullname::String
end
function PressurePiecewiseEOS(equations, transition_values)
    fullname = join([n.fullname for n in equations], " & ")
    PressurePiecewiseEOS(equations, transition_values, fullname)
end

""" Choose the appropriate individual EOS from a `PiecewiseEOS`

    * `eos`: `PiecewiseEOS` to evaluate
    * `x`: Value to evaluate at (mass or pressure, as appropriate)

    Returns a `SingleEOS`. """
function get_layer_eos(eos::PiecewiseEOS, x)
    # get layer number
    layer_edges = eos.transition_values
    for (layer_num, layer_end) in enumerate(layer_edges[2:end])
        if x < 0
            # inner layer boundary
            return eos.equations[1]::SingleEOS
        elseif x <= layer_end
            return eos.equations[layer_num]::SingleEOS
        elseif x > layer_edges[end]
            # outer layer boundary
            return eos.equations[end]::SingleEOS
        end
    end
end

import Base.length
length(eos::PiecewiseEOS) = length(eos.equations)
length(::SingleEOS) = 1

# Individual EOSes
#------------------------------------------------------------------------------

# Analytic EOSes
mgsio3_func(P) = 4100. + 0.00161*(P^0.541)
fe_func(P) = 8300. + 0.00349*(P^0.528)
h2o_func(P) = 1460. + 0.00311*(P^0.513)
graphite_func(P) = 2250. + 0.00350*(P^0.514)
sic_func(P) = 3220. + 0.00172*(P^0.537)

const mgsio3 = SimpleEOS(NoTemp, mgsio3_func, "MgSiO3")
const fe = SimpleEOS(NoTemp, fe_func, "Fe")
const h2o = SimpleEOS(NoTemp, h2o_func, "H2O")
const graphite = SimpleEOS(NoTemp, graphite_func, "Graphite")
const sic = SimpleEOS(NoTemp, sic_func, "SiC")

""" The Birch-Murnaghan EOS function, in SI units

    * `ρ`: Density
    * `ρ₀`: Reference density
    * `K₀`: Bulk modulus
    * `dK₀`: First derivative of the bulk modulus (dimensionless)

    Optional parameters
    -------------------
    * `d2K₀`: Second derivative of the bulk modulus, for 4th order BME
    """; :BME
function BME(ρ, ρ₀, K₀, dK₀)
    eta = ρ/ρ₀

    P3 = (3/2*K₀*(eta^(7/3) - eta^(5/3))
          * (1 + 3/4*(dK₀ - 4)*(eta^(2/3) - 1)))

    P3
end
function BME(ρ, ρ₀, K₀, dK₀, d2K₀)
    eta = ρ/ρ₀

    P3 = BME(ρ, ρ₀, K₀, dK₀)

    P4 = P3 + (3/2*K₀*(eta^(7/3) - eta^(5/3))
               * 3/8*(eta^(2/3) - 1)^2
               * (K₀*d2K₀ + dK₀*(dK₀ - 7) + 143/9))

    P4
end

""" The Vinet EOS function in SI units

    * `ρ`: Density
    * `ρ₀`: Reference density
    * `K₀`: Bulk modulus
    * `dK₀`: First derivative of the bulk modulus (dimensionless) """
function Vinet(ρ, ρ₀, K0, dK₀)
    eta = ρ/ρ₀

    P = (3K0*eta^(2/3) * (1 - eta^(-1/3))
         * exp(3/2*(dK₀ - 1)*(1 - eta^(-1/3))))

    P
end

""" Thomas-Fermi-Dirac EOS with energy correction in SI units.

    * `P`: Pressure
    * `Z`: Atomic number
    * `A`: Atomic mass

    To calculate the TFD for multi-atom molecules, supply `Z` and `A` as
    vectors. You may also supply `n`, the numbers of each atom, as a final
    parameter. If you do not supply `n` it will be assumed to be [1, 1, ...]
    (that is, equal numbers of all atoms). The length of `n` must match `Z` and
    `A`. """
function TFD{N<:Integer}(P, Z::Vector{N}, A::Vector, n::Vector)

    # inputs should be the same size
    @assert length(Z) == length(A) == length(n)

    # P is in Pa but we want it in dyn/cm**2: 1 Pa = 10 dyn/cm**2
    P *= 10

    # constants
    g = [0          0          0         0         0;
         0          0          0         0         0;
         1.512E-2   8.955E-2   1.090E-1  5.089     -5.980;
         2.181E-3   -4.015E-1  1.698     -9.566    9.873;
         -3.328E-4  5.167E-1   -2.369    1.349E1   -1.427E1;
         -1.384E-2  -6.520E-1  3.529     -2.095E1  2.264E1]

    # pre-calculations
    ζ   = (P / 9.524E13)^(1/5) .* Z.^(-2/3)
    ε   = (3 ./(32π^2 .* Z.^2)).^(1/3)
    ϕ   = (3^(1/3))/20 + ε./(4 .*(3 .^(1/3)))
    α   = 1./(1.941E-2 - ε.^(1/2).*6.277E-2 + ε.*1.076)
    x₀0 = (8.884E-3 + (ε.^(1/2)).*4.998E-1
                        + ε.*5.2604E-1).^(-1)
    β₀  = x₀0.*ϕ - 1
    β₁  = β₀.*α + ((1 + β₀)./ϕ)

    function β_(n::Integer)
        # sub-function for β remainder of components
        n += 1 # adjust n from 2-5 to 3-6
        bn = 1./((g[n, 1] + g[n, 2].*(ε.^(1/2)) + g[n, 3].*ε
                  + g[n, 4].*(ε.^(3/2)) + g[n, 5].*(ε.^2)).^(n-1))
    end

    β₂ = β_(2)
    β₃ = β_(3)
    β₄ = β_(4)
    β₅ = β_(5)

    βζ = β₀ + β₁.*ζ + β₂.*ζ.^2 + β₃.*ζ.^3 + β₄.*ζ.^4 + β₅.*ζ.^5

    x₀ = 1./(ζ + ϕ) .* (1 + exp(-α.*ζ).*βζ)

    num = sum(n.*A)
    denom = sum(n.*x₀.^3 ./ Z)
    ρ = num/denom * 3.866

    # rho is in g/cm3 but we want it in kg/m3: 1 g/cm3 = 1000 kg/m3
    1000ρ
end
function TFD{N<:Integer}(P, Z::Vector{N}, A::Vector)
    n = ones(A)
    TFD(P, Z, A, n)
end
function TFD(P, Z::Integer, A)
    TFD(P, [Z], [A])
end

# Interpolation, storage and retrieval
# -----------------------------------------------------------------------------
function Base.write(eos::EOS, pressures::Vector{Float64})
    densities = map(eos, pressures)
    fullname = eos.fullname
    pathname = joinpath("data", "$(fullname).eos")
    open(pathname, "w") do f
        write(f, "# Pressure\tDensity\n")
        writedlm(f, hcat(pressures, densities))
    end
end

"Put a number of interpolated EOSes into data files for later use"
function write_eoses_to_files()
    # Seager's versions of the BME
    fe_eps_seager_func(ρ) = Vinet(ρ, 8300., 156.2, 6.08) * 1e9
    h2o_VII_seager_func(ρ) = BME(ρ, 1460., 23.7, 4.15) * 1e9
    mgsio3_pv_seager_func(ρ) = BME(ρ, 4100., 247., 3.97, -0.016) * 1e9
    fe_seager_low = InvPressureEOS(fe_eps_seager_func, 1e3, 1e14,
                                "Fe (Vinet) (Seager 2007)")
    h2o_seager_low = InvPressureEOS(h2o_VII_seager_func, 1e3, 1e8,
                                 "H2O (BME3) (Seager 2007)")
    h2o_seager_dft = load_interpolated_eos("data/tabulated/H2O (DFT).eos", linear=true)
    mgsio3_seager_low = InvPressureEOS(mgsio3_pv_seager_func, 1e3, 5e4,
                                     "MgSiO3 (BME4) (Seager 2007)")

    fe_tfd_func(P) = TFD(P, 26, 55.845)
    h2o_tfd_func(P) = TFD(P, [1, 8], [1.00794, 15.9994], [2., 1.])
    mgsio3_tfd_func(P) = TFD(P, [12, 14, 8],
                                   [24.305, 28.0855, 15.9994], [1., 1., 3.])

    fe_tfd = SimpleEOS(notemp, fe_tfd_func, "Fe TFD")
    h2o_tfd = SimpleEOS(notemp, h2o_tfd_func, "H2O TFD")
    mgsio3_tfd = SimpleEOS(notemp, mgsio3_tfd_func, "MgSiO3 TFD")

    fe_seager = PressurePiecewiseEOS([fe_seager_low, fe_tfd],
                                     [1e4, 2.09e13, 1e20])
    h2o_seager = PressurePiecewiseEOS([h2o_seager_low,
                                       h2o_seager_dft,
                                       h2o_tfd],
                                      [1e4,
                                      44.3e9,
                                      7686e9, 1e20])
    mgsio3_seager = PressurePiecewiseEOS([mgsio3_seager_low, mgsio3_tfd],
                                         [1e4, 1.35e13, 1e20])

    for eos in [fe_seager, h2o_seager, mgsio3_seager]
        write(eos, logspace(4,20,500))
    end
end

""" Read a previously-written temperature-independent EOS

    If `linear`=`true`, the EOS is assumed to be on a linear grid.
    (However, the grid does not have to be regular.) """
function load_interpolated_eos(file::String; linear=false)
    data = readdlm(file, Float64; skipstart=1)
    P, rho = data[:, 1], data[:, 2]
    if linear
        interp_func = lininterp(P, rho)
    else
        interp_func = loginterp(P, rho)
    end
    _, filename = splitdir(file)
    name, _ = splitext(filename)

    SimpleEOS(NoTemp, interp_func, name)
end

""" Read a previously-written 2D temperature-independent EOS

    If `linear`=`true`, the EOS is assumed to be on a linear grid.
    However, the grid does not have to be regular. """
function load_2D_eos(file::String; linear::Bool=false)
    data = readdlm(file, Float64)
    P = vec(data[2:end, 1]) # dimension 1 (columns)
    T = vec(data[1, 2:end]) # dimension 2 (rows)
    rho = data[2:end, 2:end]

    if linear
        interp_func = lininterp(P, T, rho)
    else
        interp_func = loginterp(P, T, rho)
    end

    _, filename = splitdir(file)
    name, _ = splitext(filename)

    SimpleEOS(WithTemp, interp_func, name)
end

# EOS evaluation
# -----------------------------------------------------------------------------
import Base.call

# Inverting EOSes which are written as P(ρ)
function call(eos::InvPressureEOS, P::Real)
    fzero(x -> eos.equation(x) - P, eos.a, eos.b)
end

# Split piecewise EOSes appropriately
call(eos::PressurePiecewiseEOS, vs::ValueSet) = get_layer_eos(eos, pressure(vs))(vs)
call(eos::MassPiecewiseEOS, vs::ValueSet) = get_layer_eos(eos, mass(vs))(vs)
call(eos::PressurePiecewiseEOS, P::Real) = eos(ValueSet(0, 0, P))

# Calling EOSes with ValueSets
call(eos::EOS, vs::ValueSet) = eos(pressure(vs))
call(eos::EOS{WithTemp}, vs::PhysicalValues) = eos(pressure(vs), temperature(vs))
call(eos::EOS, P::Real) = eos.equation(P)
call(eos::EOS{WithTemp}, P::Real, T::Real) = eos.equation(P, T)

# Load interpolated EOS from file
const fe_seager = load_interpolated_eos("$DATADIR/Fe (Vinet) (Seager 2007) & Fe TFD.eos")
const h2o_seager = load_interpolated_eos("$DATADIR/H2O (BME3) (Seager 2007) & H2O (DFT) & H2O TFD.eos")
const h2o_seager_simple = load_interpolated_eos("$DATADIR/H2O (BME3) (Seager 2007) & H2O TFD.eos")
const mgsio3_seager = load_interpolated_eos("$DATADIR/MgSiO3 (BME4) (Seager 2007) & MgSiO3 TFD.eos")

const my_h2o_300 = load_interpolated_eos("$DATADIR/tabulated/h2o-300K.dat")
const my_h2o_500 = load_interpolated_eos("$DATADIR/tabulated/h2o-500K.dat")
const my_h2o_800 = load_interpolated_eos("$DATADIR/tabulated/h2o-800K.dat")
const my_h2o_1200 = load_interpolated_eos("$DATADIR/tabulated/h2o-1200K.dat")

const my_h2o_full = load_2D_eos("$DATADIR/tabulated/my_h2o_100x100.dat")

# Phase boundaries
# ------------------------------------------------------------------------------

# the 'let' bindings in this section allow us to avoid polluting the module
# scope with intermediate values, like below: we won't see 'keys', 'shortkeys',
# and so on once they've ben used to generate the phase dictionary

"A dictionary mapping phase names to short codes, like 'L' for liquid"
const phase_mappings = let
    keys = ["liquid", "ice I", "ice II", "ice III",
            "ice V", "ice VI", "ice VII", "ice VIII", "ice X"]
    shortkeys = split("L I II III V VI VII VIII X")

    phasemap = Dict(zip(keys, shortkeys))
    shortmap = Dict(zip(shortkeys, shortkeys))

    phase_mappings = merge!(phasemap, shortmap)
end

"Phase boundary parameter table from Dunaeva et al"
const dunaeva_phase_boundary_table = let 
    readdlm("$DATADIR/tabulated/Dunaeva-phase-boundaries.dat")
end

"Describes the T/P extent and shape parameters of a phase boundary"
immutable PhaseBoundaryPars
    phase1::ASCIIString
    phase2::ASCIIString
    Tmin::Float64
    Tmax::Float64
    Pmin::Float64
    Pmax::Float64
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
end
"Generate `PhaseBoundaryPars` for two given phases"
function PhaseBoundaryPars(phase1::String, phase2::String)
    table = dunaeva_phase_boundary_table

    match11 = (table[:, 1] .== phase1)
    match22 = (table[:, 2] .== phase2)
    match21 = (table[:, 2] .== phase1)
    match12 = (table[:, 1] .== phase2)
    matches_normal = match11 & match22
    matches_reversed = match12 & match21
    matches = matches_normal | matches_reversed

    row = table[matches, :]
    row[5:6] *= 1e5 # convert from bar to Pa

    PhaseBoundaryPars(row...)
end

"Holds details about the boundary between two phases"
abstract PhaseBoundary

"A phase boundary following the formulation of Dunaeva et al"
immutable DunaevaPhaseBoundary <: PhaseBoundary
    P::Vector{Float64}
    T::Vector{Float64}
    pars::PhaseBoundaryPars
end
PhaseBoundary(P::AbstractVector, T, pars) = DunaevaPhaseBoundary(collect(P), T, pars)

"Minimum pressure of the phase boundary (padded to avoid problems at P=0)"; :Pmin
Pmin(pbp::PhaseBoundaryPars) = pbp.Pmin + 1e-9
Pmin(dpb::DunaevaPhaseBoundary) = Pmin(dpb.pars)
"Maximum pressure of the phase boundary"; :Pmax
Pmax(pbp::PhaseBoundaryPars) = pbp.Pmax
Pmax(dpb::DunaevaPhaseBoundary) = Pmax(dpb.pars)

"A phase boundary that doesn't follow the Dunaeva formulation"
immutable OtherPhaseBoundary <: PhaseBoundary
    P::Vector{Float64}
    T::Vector{Float64}
end

"Get the boundary between two phases"; :PhaseBoundary
function PhaseBoundary(phase1::String, phase2::String)
    pars = PhaseBoundaryPars(phase1, phase2)
    P = linspace(Pmin(pars), Pmax(pars))
    f = P -> phase_boundary_temp(P, pars)
    T = map(f, P)

    PhaseBoundary(P, T, pars)
end
function PhaseBoundary(phase1::Symbol, phase2::Symbol)
    p1 = phase_mappings[string(phase1)]
    p2 = phase_mappings[string(phase2)]
    PhaseBoundary(p1, p2)
end
function PhaseBoundary(phase1::Symbol, phase2::String)
    p1 = phase_mappings[string(phase1)]
    PhaseBoundary(p1, phase2)
end
function PhaseBoundary(phase1::String, phase2::Symbol)
    p2 = phase_mappings[string(phase2)]
    PhaseBoundary(phase1, p2)
end

"Calculate the temperature along a phase boundary"
function phase_boundary_temp(P, pars::PhaseBoundaryPars)
    # P in Pa; this function is in bar
    Pa_to_bar(P) = P / 100000
    
    P = Pa_to_bar(P)
    T = pars.a + pars.b*P + pars.c*log(P) + pars.d/P + pars.e*sqrt(P)
end

"All relevant phase boundaries"
const phase_boundaries = let
    dunaeva_boundaries = maprows(dunaeva_phase_boundary_table) do row
        p1, p2 = row[1:2]
        PhaseBoundary(p1, p2)
    end

    iapws_boundary = let
        table = readdlm("$DATADIR/tabulated/iapws-phase-boundary.dat")
        P = table[:, 1] * 1e6 # MPa -> Pa
        T = table[:, 2]
        OtherPhaseBoundary(P, T)
    end

    vcat(dunaeva_boundaries, iapws_boundary)
end
