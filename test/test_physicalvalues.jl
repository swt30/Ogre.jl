using Base.Test
import Ogre


@testset "Value sets" begin
    vs_notemp = Ogre.ValueSet(1,2,3)
    vs_full = Ogre.ValueSet(1,2,3,4)
    @test_throws MethodError Ogre.ValueSet(1,2)

    @testset "Type hierarchy" begin
        tempdep = Ogre.ValueSet{Ogre.WithTemp}
        notemp = Ogre.ValueSet{Ogre.NoTemp}
        @test !isa(vs_notemp, tempdep)
        @test isa(vs_notemp, notemp)
        @test isa(vs_notemp, Ogre.ValueSet)
        @test isa(vs_full, tempdep)
        @test !isa(vs_full, notemp)
        @test isa(vs_full, Ogre.ValueSet)
    end

    @testset "Get values from them" begin
        @test Ogre.mass(vs_notemp) == 1
        @test Ogre.nonmass(vs_notemp) == [2,3]
        @test Ogre.mass(vs_full) == 1
        @test Ogre.nonmass(vs_full) == [2,3,4]
        @test Ogre.radius(vs_notemp) == 2
        @test Ogre.pressure(vs_full) == 3
        @test_throws MethodError Ogre.temperature(vs_notemp)
        @test Ogre.temperature(vs_full) == 4
    end

    @testset "Zero values or less are considered 'unphysical'" begin
        @test !Ogre.isphysical(Ogre.zero(Ogre.ValueSet{Ogre.NoTemp}))
        @test !Ogre.isphysical(Ogre.zero(Ogre.ValueSet{Ogre.WithTemp}))
        @test !Ogre.isphysical(Ogre.ValueSet(0, 0, -1))
        @test !Ogre.isphysical(Ogre.ValueSet(0, 0, 0, -1))
    end
end
