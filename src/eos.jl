# EOS.JL
# Equation of state handling

# Required packages
using Dierckx, Roots

# EOS types
#------------------------------------------------------------------------------

abstract EOS <: Equation
abstract SingleEOS <: EOS
abstract PiecewiseEOS <: EOS
abstract SimpleEOS <: SingleEOS

@doc """
    Type for simple equations of state (one function).

    * `equation`: Function ρ=f(P)
    * `fullname`: Name of the EOS (for printing and plots)
    """ ->
immutable PressureEOS <: SimpleEOS
    equation::Function
    fullname::String
end

@doc """
    Type for simple equations of state (one function).

    * `equation`: Function ρ=f(P, T)
    * `fullname`: Name of the EOS (for printing and plots)
    """ ->
immutable PressureTempEOS <: SimpleEOS
    equation::Function
    fullname::String
end

@doc """
    Type for equations of state which must invert some function over a
    density range.

    * `equation`: Function P=f(ρ)
    * `a`, `b`: Density range to invert over
    * `fullname`: Name of the EOS (for printing and plots)
    """ ->
immutable InvPressureEOS{T<:Real} <: SingleEOS
    equation::Function
    a::T
    b::T
    fullname::String
end

@doc """
    Type for equations of state which are piecewise in the mass coordinate.

    * `equations`: Vector of `SingleEOS` for each piece.
    * `transition_values`: Vector defining the edge of each piece. Evaluating
      past the ends of this vector will just evaluate the nearest EOS.
    * `fullname`: Name of the EOS (for printing and plotting)
    """ ->
immutable MassPiecewiseEOS{T<:Real, E<:SingleEOS} <: PiecewiseEOS
    equations::Vector{E}
    transition_values::Vector{T}
    fullname::String
end
function MassPiecewiseEOS{T<:Real, E<:SingleEOS}(equations::Vector{E},
    transition_values::Vector{T})

    fullname = join([n.fullname for n in equations], " & ")
    MassPiecewiseEOS(equations, transition_values, fullname)
end
function MassPiecewiseEOS{T<:Real, E<:SingleEOS}(equations::Vector{E},
    M::Real, mass_fractions::Vector{T})

    layer_edges = [0, cumsum(mass_fractions)] .* M
    MassPiecewiseEOS(equations, layer_edges)
end

@doc """
    Type for equations of state which are piecewise in the pressure coordinate.

    * `equations`: Vector of `SingleEOS` for each piece.
    * `transition_values`: Vector defining the edge of each piece. Evaluating
      past the ends of this vector will just evaluate the nearest EOS.
    * `fullname`: Name of the EOS (for printing and plotting)
    """ ->
immutable PressurePiecewiseEOS{T<:Real, E<:SingleEOS} <: PiecewiseEOS
    equations::Vector{E}
    transition_values::Vector{T}
    fullname::String
end
function PressurePiecewiseEOS{T<:Real, E<:SingleEOS}(equations::Vector{E},
    transition_values::Vector{T})

    fullname = join([n.fullname for n in equations], " & ")
    PressurePiecewiseEOS(equations, transition_values, fullname)
end

@doc """
    Get the appropriate individual EOS from a `PiecewiseEOS`, whether
    mass-piecewise or pressure-piecewise.

    * `eos`: `PiecewiseEOS` to evaluate
    * `x`: Value to evaluate at

    Returns a `SingleEOS`.
    """ ->
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

@doc """Get the number of equations in an `EOS`""" ->
n_eqs(eos::PiecewiseEOS) = length(eos.equations)
n_eqs(::SingleEOS) = 1

# Individual EOSes
#------------------------------------------------------------------------------

# Analytic EOS

mgsio3_func(P::Real) = 4100. + 0.00161*(P^0.541)
fe_func(P::Real) = 8300. + 0.00349*(P^0.528)
h2o_func(P::Real) = 1460. + 0.00311*(P^0.513)
graphite_func(P::Real) = 2250. + 0.00350*(P^0.514)
sic_func(P::Real) = 3220. + 0.00172*(P^0.537)

mgsio3 = PressureEOS(mgsio3_func, "MgSiO3")
fe = PressureEOS(fe_func, "Fe")
h2o = PressureEOS(h2o_func, "H2O")
graphite = PressureEOS(graphite_func, "Graphite")
sic = PressureEOS(sic_func, "SiC")

@doc """
    The Birch-Murnaghan EOS function, in SI units

    * `rho`: Density
    * `rho0`: Reference density
    * `K0`: Bulk modulus
    * `dK0`: First derivative of the bulk modulus (dimensionless)

    Optional parameters
    -------------------
    * `d2K0`: Second derivative of the bulk modulus for 4th order BME
    """ ->
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
BME(args...) = BME(promote(args)...)

@doc """
    The Vinet EOS function in SI units

    * `rho`: Density
    * `rho0`: Reference density
    * `K0`: Bulk modulus
    * `dK0`: First derivative of the bulk modulus (dimensionless)
    """ ->
function Vinet{T<:Real}(rho::T, rho0::T, K0::T, dK0::T)
    eta = rho/rho0

    P = (3K0*eta^(2/3) * (1 - eta^(-1/3))
         * exp(3/2*(dK0 - 1)*(1 - eta^(-1/3))))

    P
end
Vinet(args...) = Vinet(promote(args)...)

@doc """
    Thomas-Fermi-Dirac EOS with energy correction in SI units.

    * `P`: Pressure
    * `Z`: Atomic number
    * `A`: Atomic mass

    To calculate the TFD for multi-atom molecules, supply `Z` and `A` as
    vectors. You may also supply `n`, the numbers of each atom, as a final
    parameter. If you do not supply `n` it will be assumed to be [1, 1, ...]
    (equal numbers of all atoms). The length of `n` must match `Z` and `A`.
    """ ->
function TFD{T<:Real, N<:Integer}(P::T,
    Z::Vector{N}, A::Vector{T}, n::Vector{T})

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
    1000ρ::T
end
function TFD{T<:Real, N<:Integer}(P::T, Z::Vector{N}, A::Vector{T})
    n = ones(A)
    TFD(P, Z, A, n)
end
function TFD{T<:Real}(P::T, Z::Integer, A::T)
    TFD(P, [Z], [A])
end
TFD(args...) = TFD(promote(args)...)

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

@doc """Put a number of interpolated EOSes into data files for later use""" ->
function write_eoses_to_files()
    # Seager's versions of the BME
    fe_eps_seager_func(rho::Real) = Vinet(rho, 8300., 156.2, 6.08) * 1e9
    h2o_VII_seager_func(rho::Real) = BME(rho, 1460., 23.7, 4.15) * 1e9
    mgsio3_pv_seager_func(rho::Real) = BME(rho, 4100., 247., 3.97, -0.016) * 1e9
    fe_seager_low = InvertedEOS(fe_eps_seager_func, 1e3, 1e14,
                                "Fe (Vinet) (Seager 2007)")
    h2o_seager_low = InvertedEOS(h2o_VII_seager_func, 1e3, 1e8,
                                 "H2O (BME3) (Seager 2007)")
    h2o_seager_dft = load_interpolated_eos("data/tabulated/H2O (DFT).eos", linear=true)
    mgsio3_seager_low = InvertedEOS(mgsio3_pv_seager_func, 1e3, 5e4,
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
                                      [1e4,
                                      44.3e9,
                                      7686e9, 1e20])
    mgsio3_seager = PressurePiecewiseEOS([mgsio3_seager_low, mgsio3_tfd],
                                         [1e4, 1.35e13, 1e20])

    for eos in [fe_seager, h2o_seager, mgsio3_seager]
        write(eos, logspace(4,20,500))
    end
end

@doc """
    Create an interpolating function using a log-spaced coordinate grid.
    """ ->
function loginterp{T<:Real}(x::Array{T}, y::Array{T})
    # first transform the grid to be linear
    logx = log10(x)
    # then do the interpolation as if it were linear
    lin_interp_func = lininterp(logx, y)

    interp_func(x) = lin_interp_func(log10(x))
end

@doc """Create an interpolating function in linear space""" ->
function lininterp{T<:Real}(x::Array{T}, y::Array{T})

    spline = Spline1D(x, y, k=2)

    interp_func(x) = evaluate(spline, convert(Float64, x))
end

@doc """
    Read a previous-written EOS into a `SimpleEOS`

    If `linear`=`true`, the EOS is assumed to be on a linear grid.
    However, the grid does not have to be regular.
    """ ->
function load_interpolated_eos(file::String; linear=false)
    data = readdlm(file, Float64; skipstart=1)
    P, rho = data[:, 1], data[:, 2]
    if linear
        interp_func = lininterp(P, rho)
    else
        interp_func = loginterp(P, rho)
    end
    directory, filename = splitdir(file)
    name, ext = splitext(filename)

    PressureEOS(interp_func, name)
end

# EOS evaluation
# -----------------------------------------------------------------------------
import Base.call

# Inverting EOSes which are written as P(ρ)
function call(eos::InvPressureEOS, P::Real)
    fzero(x -> eos.equation(x) - P, eos.a, eos.b)
end

# Split piecewise EOSes appropriately
call(eos::PressurePiecewiseEOS, vs::ValueSet) = get_layer_eos(eos, vs.P)(vs)
call(eos::MassPiecewiseEOS, vs::ValueSet) = get_layer_eos(eos, vs.m)(vs)
call(eos::PressurePiecewiseEOS, P::Real) = eos(ValueSet(0, 0, P))

# Calling EOSes with ValueSets
call(eos::SingleEOS, vs::MassRadiusPressure) = eos(vs.P)
call(eos::SingleEOS, vs::PhysicalValues) = eos(vs.P, vs.T)
call(eos::PressureEOS, P::Real) = eos.equation(P)
call(eos::PressureTempEOS, P::Real, T::Real) = eos.equation(P, T)

# Load interpolated EOS from file
fe_seager = load_interpolated_eos("$DATADIR/Fe (Vinet) (Seager 2007) & Fe TFD.eos")
h2o_seager = load_interpolated_eos("$DATADIR/H2O (BME3) (Seager 2007) & H2O (DFT) & H2O TFD.eos")
h2o_seager_simple = load_interpolated_eos("$DATADIR/H2O (BME3) (Seager 2007) & H2O TFD.eos")
mgsio3_seager = load_interpolated_eos("$DATADIR/MgSiO3 (BME4) (Seager 2007) & MgSiO3 TFD.eos")

my_h2o_300 = load_interpolated_eos("$DATADIR/tabulated/h2o-300K.dat")
my_h2o_500 = load_interpolated_eos("$DATADIR/tabulated/h2o-500K.dat")
my_h2o_800 = load_interpolated_eos("$DATADIR/tabulated/h2o-800K.dat")
my_h2o_1200 = load_interpolated_eos("$DATADIR/tabulated/h2o-1200K.dat")

