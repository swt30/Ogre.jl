# Tests on utility functions

using FactCheck
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


facts("Utility function tests") do
    res = test_util_resources

    context("Copy-modification of types") do
        t1 = res.my_immutable_type("hello", 12)
        t2 = Ogre.cpmod(t1, field1="world")
        @fact t2.field1 --> "world"
        @fact t2.field2 --> 12

        t3 = res.my_mutable_type("foo", 13)
        t4 = Ogre.cpmod(t3, field2=14)
        @fact t4.field1 --> "foo"
        @fact t4.field2 --> 14
    end

    context("Mapping rows of a function") do
        M = [1 2 3;
             4 5 6;
             7 8 9]
        @fact vec(Ogre.maprows(sum, M)) --> [6, 15, 24]
        @fact vec(Ogre.maprows(mean, M)) --> [2., 5., 8.]
    end

    context("Testing for NaN values") do
        @fact Ogre.hasnan([1, 2, 3, 4]) --> false
        @fact Ogre.hasnan([1, NaN, 3, NaN]) --> true
    end
end
