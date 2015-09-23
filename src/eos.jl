# Equations of state

using WaterData
using Iterators: chain

# The bulk of the EOSes are provided in the WaterData package.
# We define one more EOS which is mass-piecewise (for modelling layered planets)
""" Equation of state which is piecewise in the mass coordinate

    * `eoses`: Vector of EOSes for each piece.
    * `edges`: Vector defining the edge of each piece. Evaluating
      past the ends of this vector will just evaluate the nearest EOS. """
immutable MassPiecewiseEOS{E<:EOS, T<:Real} <: WaterData.PiecewiseEOS
    eoses::Vector{E}
    edges::Vector{T}
end
function MassPiecewiseEOS(eoses, M, mass_fractions)
    layer_edges = vcat([0], cumsum(mass_fractions)) .* M
    MassPiecewiseEOS(eoses, layer_edges)
end


# EOS evaluation
# (we loop through a bunch of EOSes to avoid ambiguity clashes)

# Split piecewise EOSes appropriately
for T in (MassRadiusPressure, PhysicalValues)
    Base.call(eos::MassPiecewiseEOS, vs::T) = WaterData.get_single_eos(eos, mass(vs))(vs)
end

# Calling EOSes with ValueSets
let specific_eoses = (TFD, BME3, BME4, Vinet, PolytropicEOS),
    modifier_eoses = (BoundedEOS, WaterData.InverseFunctionalEOS, PressurePiecewiseEOS),
    special_eoses = (WaterData.OutOfDomainEOS, ),
    general_eoses = (LineEOS, GridEOS, EOS)

    for T in chain(specific_eoses, special_eoses, modifier_eoses, general_eoses)
        Base.call(eos::T, vs::MassRadiusPressure) = eos(pressure(vs))
        Base.call(eos::T, vs::PhysicalValues) = eos(pressure(vs), temperature(vs))
    end
end
