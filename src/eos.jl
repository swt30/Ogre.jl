module eos
export EOS, SimpleEOS, PiecewiseEOS, mgsio3, fe, h2o, callfunc
using ogre.common, Grid

abstract EOS <: Equation

# Equations of state

type SimpleEOS <: EOS
    equation::Function
    name::String
end

type PiecewiseEOS <: EOS
    equations::Vector{SimpleEOS}
    transition_masses::Vector{Float64}
end

function get_layer_eos(eos::PiecewiseEOS, m::Real)
    # get layer number
    layer_edges = eos.transition_masses
    for (layer_num, layer_end) in enumerate(layer_edges[2:end])
        if m < 0
            # inner layer boundary
            println("m became smaller than 0")
            return eos.equations[1]::SimpleEOS
        elseif m <= layer_end
            return eos.equations[layer_num]::SimpleEOS
        elseif m > layer_edges[end]
            # outer layer boundary
            return eos.equations[end]::SimpleEOS
        end
    end
end

n_eqs(eos::PiecewiseEOS) = length(eos.equations)
n_eqs(eos::SimpleEOS) = 1

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

# my tabular H2O EOS

f = readdlm("$(datadir())/h2o.dat")

x = vec(f[:, 1])
a = x[1]
b = x[end]
y = vec(f[:, 2])

logrange = linrange(log10(a), log10(b), length(x))
yi = CoordInterpGrid(logrange, y, BCnan, InterpQuadratic)
my_h2o_func(P) = yi[log10(P)]

my_h2o = SimpleEOS(my_h2o_func)

# density retrieval functions for EOS in particular
import ogre.common.callfunc

function callfunc(eos::SimpleEOS, vs::ValueSet)
    callfunc(eos, vs.P)
end

function callfunc(eos::PiecewiseEOS, vs::ValueSet)
    single_eos = get_layer_eos(eos, vs.m)
    callfunc(single_eos, vs.P)
end

end
