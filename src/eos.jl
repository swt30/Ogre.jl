module Eos
using Ogre: Common
using Grid, Roots
# Exported types
export EOS, SimpleEOS, MassPiecewiseEOS, InvertedEOS
# Exported functions
export callfunc

# Equations of state

abstract EOS <: Equation
abstract SingleEOS <: EOS
abstract PiecewiseEOS <: EOS

immutable SimpleEOS <: SingleEOS
    equation::Function
    fullname::String
end

immutable InvertedEOS{T<:Real} <: SingleEOS
    equation::Function
    a::T
    b::T
    fullname::String
end

immutable MassPiecewiseEOS{T<:Real, E<:SingleEOS} <: PiecewiseEOS
    equations::Vector{E}
    transition_values::Vector{T}
    fullname::String
end

immutable PressurePiecewiseEOS{T<:Real, E<:SingleEOS} <: PiecewiseEOS
    equations::Vector{E}
    transition_values::Vector{T}
    fullname::String
end

# constructors for the two piecewise EOSes
function PressurePiecewiseEOS{T<:Real, E<:SingleEOS}(equations::Vector{E},
                                             transition_values::Vector{T})
    fullname = join([n.fullname for n in equations], " & ")
    PressurePiecewiseEOS(equations, transition_values, fullname)
end

function MassPiecewiseEOS{T<:Real, E<:SingleEOS}(equations::Vector{E},
                                             transition_values::Vector{T})
    fullname = join([n.fullname for n in equations], " & ")
    MassPiecewiseEOS(equations, transition_values, fullname)
end

#= get the appropriate individual EOS from a PiecewiseEOS... you can do this
either for a mass-differentiated EOS (e.g. layers in a planet) or for a
pressure-differentiated EOS (e.g. switching over from one EOS to another at
higher pressures) -- the representation is the same both ways =#
function get_layer_eos(eos::PiecewiseEOS, x::Real)
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

n_eqs(eos::PiecewiseEOS) = length(eos.equations)
n_eqs(::SimpleEOS) = 1

# polytropic EOS

mgsio3_func(P::Real) = 4100. + 0.00161*(P^0.541)
fe_func(P::Real) = 8300. + 0.00349*(P^0.528)
h2o_func(P::Real) = 1460. + 0.00311*(P^0.513)
graphite_func(P::Real) = 2250. + 0.00350*(P^0.514)
sic_func(P::Real) = 3220. + 0.00172*(P^0.537)

mgsio3 = SimpleEOS(mgsio3_func, "MgSiO3")
fe = SimpleEOS(fe_func, "Fe")
h2o = SimpleEOS(h2o_func, "H2O")
graphite = SimpleEOS(graphite_func, "Graphite")
sic = SimpleEOS(sic_func, "SiC")

# general forms of the BME and Vinet EOS

function BME{T<:Real}(rho::T, rho0::T, K0::T, dK0::T)
    eta = rho/rho0

    P3 = (3/2*K0*(eta^(7/3) - eta^(5/3))
          * (1 + 3/4*(dK0 - 4)*(eta^(2/3) - 1)))

    P3
end

function BME{T<:Real}(rho::T, rho0::T, K0::T, dK0::T, d2K0::T)
    eta = rho/rho0

    P3 = BME(rho, rho0, K0, dK0)

    P4 = P3 + (3/2*K0*(eta^(7/3) - eta^(5/3))
               * 3/8*(eta^(2/3) - 1)^2
               * (K0*d2K0 + dK0*(dK0 - 7) + 143/9))

    P4
end

function Vinet{T<:Real}(rho::T, rho0::T, K0::T, dK0::T)
    eta = rho/rho0

    P = (3K0*eta^(2/3) * (1 - eta^(-1/3))
         * exp(3/2*(dK0 - 1)*(1 - eta^(-1/3))))

    P
end

# Here is the Python TFD

# using PyCall
# eosfuncs = pyimport("eos.funcs")
# pyTFD = eosfuncs[:TFD]

# function TFD{T<:Real, N<:Integer}(P::T, Z::Vector{N}, A::Vector{T},
#                                   n::Vector{T}; python=true)
#     P *= 10
#     rho::Float64 = pyTFD(P, Z, A, n)[1]
#     rho * 1000
# end

# Here's the real, working TFD

function TFD{T<:Real, N<:Integer}(P::T, Z::Vector{N}, A::Vector{T},
                                  n::Vector{T})
    # inputs should be the same size
    @assert length(Z) == length(A)
    @assert length(n) == length(A)

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
        """sub-function for β remainder of components"""
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
    1000ρ::T
end

function TFD{T<:Real, N<:Integer}(P::T, Z::Vector{N}, A::Vector{T})
    n = ones(A)
    TFD(P, Z, A, n)
end

function TFD{T<:Real}(P::T, Z::Integer, A::T)
    TFD(P, [Z], [A])
end

# EOS interpolation, storage and retrieval

import Base.write
function write(eos::EOS, pressures::Vector{Float64})
    densities = map(pressures) do P
        callfunc(eos, P)
    end
    fullname = eos.fullname
    pathname = joinpath("data", "$(fullname).eos")
    open(pathname, "w") do f
        write(f, "# Pressure\tDensity\n")
        writedlm(f, hcat(pressures, densities))
    end
end

function write_eoses_to_files()
    # Seager's versions of the BME
    fe_eps_seager_func(rho::Real) = Vinet(rho, 8300., 156.2, 6.08) * 1e9
    h2o_VII_seager_func(rho::Real) = BME(rho, 1460., 23.7, 4.15) * 1e9
    mgsio3_pv_seager_func3(rho::Real) = BME(rho, 4100., 247., 3.97) * 1e9
    mgsio3_pv_seager_func4(rho::Real) = BME(rho, 4100., 247., 3.97, -0.016) * 1e9
    fe_seager_low = InvertedEOS(fe_eps_seager_func, 1e3, 1e14,
                                "Fe (Vinet) (Seager 2007)")
    h2o_seager_low = InvertedEOS(h2o_VII_seager_func, 1e3, 1e8,
                                 "H2O (BME3) (Seager 2007)")
    h2o_seager_dft = load_interpolated_eos("data/tabulated/H2O (DFT).eos")
    # mgsio3_seager3_low = InvertedEOS(mgsio3_pv_seager_func3, 1e3, 1e8,
    #                                  "MgSiO3 (BME3) (Seager 2007)")
    mgsio3_seager4_low = InvertedEOS(mgsio3_pv_seager_func4, 1e3, 5e4,
                                     "MgSiO3 (BME4) (Seager 2007)")

    fe_tfd_func(P::Real) = TFD(P, 26, 55.845)
    h2o_tfd_func(P::Real) = TFD(P, [1, 8], [1.00794, 15.9994], [2., 1.])
    mgsio3_tfd_func(P::Real) = TFD(P, [12, 14, 8],
                                   [24.305, 28.0855, 15.9994], [1., 1., 3.])

    fe_tfd = SimpleEOS(fe_tfd_func, "Fe TFD")
    h2o_tfd = SimpleEOS(h2o_tfd_func, "H2O TFD")
    mgsio3_tfd = SimpleEOS(mgsio3_tfd_func, "MgSiO3 TFD")

    fe_seager = PressurePiecewiseEOS([fe_seager_low, fe_tfd],
                                     [1e4, 2.09e13, 1e20])
    h2o_seager = PressurePiecewiseEOS([h2o_seager_low,
                                       h2o_seager_dft,
                                       h2o_tfd],
                                      [1e4, 44.3e9, 7686e9, 1e20])
    mgsio3_seager = PressurePiecewiseEOS([mgsio3_seager4_low, mgsio3_tfd],
                                         [1e4, 1.35e13, 1e20])

    for eos in [fe_seager, h2o_seager, mgsio3_seager]
        write(eos, logspace(4,20,500))
    end
end

function loginterp{T<:Real}(x::Array{T}, y::Array{T})
    min_x, max_x = extrema(x)
    logrange = linrange(log10(min_x), log10(max_x), length(x))
    yi = CoordInterpGrid(logrange, y, BCnan, InterpQuadratic)

    interp_func(P) = yi[log10(P)]
end

function lininterp{T<:Real}(x::Array{T}, y::Array{T})
    min_x, max_x = extrema(x)
    range = linrange(min_x, max_x, length(x))
    yi = CoordInterpGrid(range, y, BCnan, InterpQuadratic)

    interp_func(P) = yi[P]
end

function load_interpolated_eos(file::String)
    data = readdlm(file, Float64; skipstart=1)
    P, rho = data[:, 1], data[:, 2]
    interp_func = loginterp(P, rho)
    directory, filename = splitdir(file)
    name, ext = splitext(filename)

    SimpleEOS(interp_func, name)
end

# density retrieval functions for EOS in particular
import Ogre.Common.callfunc

function callfunc(eos::InvertedEOS, P::Real)
    fzero(x -> eos.equation(x) - P, eos.a, eos.b)
end

function callfunc(eos::SingleEOS, vs::ValueSet)
    callfunc(eos, vs.P)
end

function callfunc(eos::MassPiecewiseEOS, vs::ValueSet)
    single_eos = get_layer_eos(eos, vs.m)
    callfunc(single_eos, vs.P)
end

function callfunc(eos::PressurePiecewiseEOS, vs::ValueSet)
    single_eos = get_layer_eos(eos, vs.P)
    callfunc(single_eos, vs.P)
end

function callfunc(eos::PressurePiecewiseEOS, P::Real)
    callfunc(eos, ValueSet(0., 0., P))
end

# Generate and export interpolated functions

my_h2o = load_interpolated_eos("$DATADIR/tabulated/h2o.dat")
fe_seager = load_interpolated_eos("$DATADIR/Fe (Vinet) (Seager 2007) & Fe TFD.eos")
h2o_seager = load_interpolated_eos("$DATADIR/H2O (BME3) (Seager 2007) & H2O (DFT) & H2O TFD.eos")
mgsio3_seager = load_interpolated_eos("$DATADIR/MgSiO3 (BME4) (Seager 2007) & MgSiO3 TFD.eos")

end
