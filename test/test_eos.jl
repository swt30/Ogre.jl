import FactCheck
import Ogre


# test resources

module test_eos_resources
    import Ogre

    type SimpleEOS <: Ogre.EOS
        A::Float64
        n::Float64
        C::Float64
    end
    Ogre.@addEOSCall SimpleEOS

    type ComplicatedEOS <: Ogre.EOS
        A::Float64
        B::Float64
    end
    Ogre.@addEOSCall ComplicatedEOS

    type ConstantEOS <: Ogre.EOS
        C::Float64
    end
    Ogre.@addEOSCall ConstantEOS

    Base.call(eos::SimpleEOS, P::Number) = eos.A*P^eos.n + eos.C
    Base.call(eos::SimpleEOS, P::Number, T::Number) = eos(P)
    Base.call(eos::ComplicatedEOS, P::Number, T::Number) = eos.A*P + eos.B*T
    Base.call(eos::ConstantEOS, P::Number) = eos.C
    Base.call(eos::ConstantEOS, P::Number, T::Number) = eos(P)

    Ogre.istempdependent(::Union{SimpleEOS, ConstantEOS}) = false
    Ogre.istempdependent(::ComplicatedEOS) = true
end

facts("Basic equation of state tests") do
    res = test_eos_resources

    eos1 = res.SimpleEOS(1.,2.,3.)
    eos2 = res.ComplicatedEOS(1.,2.)
    eos3 = res.ConstantEOS(4.)

    vs = Ogre.ValueSet
    vs1 = vs(1,2,3)
    vs2 = vs(1,2,3,4)

    context("Calling with physical value sets") do
        @fact eos1(vs1) --> 1*3^2 + 3
        @fact eos1(vs2) --> 1*3^2 + 3
        @fact_throws eos2(vs1)
        @fact eos2(vs2) --> 1*3 + 2*4
        @fact eos3(vs1) --> 4
        @fact eos3(vs2) --> 4
    end

    context("Handling mass piecewise EOSes") do
        mpw = Ogre.MassPiecewiseEOS([eos1, eos2, eos3], [0, 1, 2, 3])
        mpw2 = Ogre.MassPiecewiseEOS([eos1, eos2, eos3], 3, [1/3, 1/3, 1/3])

        @fact mpw(vs(0.5, 1, 2)) --> eos1(2)
        @fact mpw2(vs(0.5, 1, 2)) --> eos1(2)
        @fact mpw(vs(1.5, 1, 2, 3)) --> eos2(2, 3)
        @fact mpw2(vs(1.5, 1, 2, 3)) --> eos2(2, 3)
        @fact mpw(vs(2.5, 1, 2, 3)) --> eos3(2, 3)
        @fact mpw2(vs(2.5, 1, 2, 3)) --> eos3(2, 3)
        @fact_throws(mpw(vs(1.5, 2, 3)))
        @fact_throws(mpw2(vs(1.5, 2, 3)))
    end
end
