module eos
export EOS, SimpleEOS, PiecewiseEOS, mgsio3, fe, h2o, call
using ogre.common, Grid

abstract EOS <: Equation

# Equations of state

type SimpleEOS <: EOS
    equation::Function
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

function mgsio3_func(P::Real)
    if P > 0
        return 4260. + 0.00127*(P^0.549)
    else
        return 4260.
    end
end

function fe_func(P::Real)
    if P > 0
        return 8300. + 0.00349*(P^0.528)
    else
        return 8300.
    end
end

mgsio3 = SimpleEOS(mgsio3_func)
fe = SimpleEOS(fe_func)

# tabular H2O EOS

f = readdlm("$(datadir())/h2o.dat")

x = vec(f[:, 1])
a = x[1]
b = x[end]
y = vec(f[:, 2])

logrange = linrange(log10(a), log10(b), length(x))
yi = CoordInterpGrid(logrange, y, BCnan, InterpQuadratic)
h2o_func(P) = yi[log10(P)]

h2o = SimpleEOS(h2o_func)

# density retrieval functions for EOS in particular
import ogre.common.call

function call(eos::SimpleEOS, vs::ValueSet)
    call(eos, vs.P)
end

function call(eos::PiecewiseEOS, vs::ValueSet)
    single_eos = get_layer_eos(eos, vs.m)
    call(single_eos, vs.P)
end

end
