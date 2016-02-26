# Equations of state

using WaterData
import JLD, Dierckx, VoronoiDelaunay, GeometricalPredicates
import WaterData: extracteos
using Iterators: chain

# The bulk of the EOSes and heat capacities are provided in the WaterData
# package. We define one more EOS which is mass-piecewise (for modelling layered
# planets)
""" Equation of state which is piecewise in the mass coordinate

    * `eoses`: Vector of EOSes for each piece.
    * `edges`: Vector defining the edge of each piece. Evaluating
      past the ends of this vector will just evaluate the nearest EOS. """
immutable MassPiecewiseEOS{E<:EOS, T<:Real} <: WaterData.PiecewiseEOS
    eoses::Vector{E}
    edges::Vector{T}
end
function MassPiecewiseEOS(eoses, M, mass_fractions)
    layer_edges = vcat([-Inf], cumsum(mass_fractions)) .* M
    layer_edges[end] = Inf
    MassPiecewiseEOS(eoses, layer_edges)
end

# EOS evaluation
# (we loop through a bunch of EOSes to avoid ambiguity clashes)

# Split piecewise EOSes appropriately
function extracteos(eos::MassPiecewiseEOS, vs::ValueSet)
    extracteos(eos, mass(vs))
end
function extracteos(eos::WaterData.PressurePiecewiseEOS, vs::ValueSet)
    extracteos(eos, pressure(vs))
end
for T in (MassRadiusPressure, PhysicalValues)
    Base.call(eos::MassPiecewiseEOS, vs::T) = extracteos(eos, vs)(vs)
end

# Calling EOSes with ValueSets
let specific_eoses = (TFD, BME3, BME4, Vinet, PolytropicEOS, WaterData.OutOfDomainEOS),
    modifier_eoses = (BoundedEOS, WaterData.InverseFunctionalEOS, PressurePiecewiseEOS),
    general_eoses = (LineEOS, GridEOS, EOS, )

    for T in chain(specific_eoses, modifier_eoses, general_eoses)
        Base.call(eos::T, vs::MassRadiusPressure) = eos(pressure(vs))
        Base.call(eos::T, vs::PhysicalValues) = eos(pressure(vs), temperature(vs))
    end
end
