using Base.Test
import Ogre


# test resources

module test_eos_resources
    import Ogre

    type SimpleEOS <: Ogre.EOS
        A::Float64
        n::Float64
        C::Float64
    end

    type ComplicatedEOS <: Ogre.EOS
        A::Float64
        B::Float64
    end

    type ConstantEOS <: Ogre.EOS
        C::Float64
    end

    for eos in (SimpleEOS, ComplicatedEOS, ConstantEOS)
        Ogre.@addEOSCall eos
    end

    (eos::SimpleEOS)(P::Number) = eos.A*P^eos.n + eos.C
    (eos::SimpleEOS)(P::Number, T::Number) = eos(P)
    (eos::ComplicatedEOS)(P::Number, T::Number) = eos.A*P + eos.B*T
    (eos::ConstantEOS)(P::Number) = eos.C
    (eos::ConstantEOS)(P::Number, T::Number) = eos(P)

    Ogre.istempdependent(::Union{SimpleEOS, ConstantEOS}) = false
    Ogre.istempdependent(::ComplicatedEOS) = true
end

@testset "Basic equation of state tests" begin
    res = test_eos_resources

    eos1 = res.SimpleEOS(1.,2.,3.)
    eos2 = res.ComplicatedEOS(1.,2.)
    eos3 = res.ConstantEOS(4.)

    vs = Ogre.ValueSet
    vs1 = vs(1,2,3)
    vs2 = vs(1,2,3,4)

    @testset "Calling with physical value sets" begin
        @test eos1(vs1) == 1*3^2 + 3
        @test eos1(vs2) == 1*3^2 + 3
        @test_throws MethodError eos2(vs1)
        @test eos2(vs2) == 1*3 + 2*4
        @test eos3(vs1) == 4
        @test eos3(vs2) == 4
    end

    @testset "Handling mass piecewise EOSes" begin
        mpw = Ogre.MassPiecewiseEOS([eos1, eos2, eos3], [0, 1, 2, 3])
        mpw2 = Ogre.MassPiecewiseEOS([eos1, eos2, eos3], 3, [1/3, 1/3, 1/3])

        @test mpw(vs(0.5, 1, 2)) == eos1(2)
        @test mpw2(vs(0.5, 1, 2)) == eos1(2)
        @test mpw(vs(1.5, 1, 2, 3)) == eos2(2, 3)
        @test mpw2(vs(1.5, 1, 2, 3)) == eos2(2, 3)
        @test mpw(vs(2.5, 1, 2, 3)) == eos3(2, 3)
        @test mpw2(vs(2.5, 1, 2, 3)) == eos3(2, 3)
        @test_throws MethodError mpw(vs(1.5, 2, 3))
        @test_throws MethodError mpw2(vs(1.5, 2, 3))
    end
end
