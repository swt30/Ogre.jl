# Tests on utility functions

if VERSION < v"0.5"
    using BaseTestNext
else
    using Base.Test
end

import Ogre

module test_util_resources
    immutable my_immutable_type
        field1
        field2
    end

    type my_mutable_type
        field1
        field2
    end
end


@testset "Utility function tests" begin
    res = test_util_resources

    @testset "Copy-modification of types" begin
        t1 = res.my_immutable_type("hello", 12)
        t2 = Ogre.cpmod(t1, field1="world")
        @test t2.field1 == "world"
        @test t2.field2 == 12

        t3 = res.my_mutable_type("foo", 13)
        t4 = Ogre.cpmod(t3, field2=14)
        @test t4.field1 == "foo"
        @test t4.field2 == 14
    end

    @testset "Mapping rows of a function" begin
        M = [1 2 3;
             4 5 6;
             7 8 9]
        @test vec(Ogre.maprows(sum, M)) == [6, 15, 24]
        @test vec(Ogre.maprows(mean, M)) == [2., 5., 8.]
    end

    @testset "Testing for NaN values" begin
        @test !Ogre.hasnan([1, 2, 3, 4])
        @test Ogre.hasnan([1, NaN, 3, NaN])
    end
end
