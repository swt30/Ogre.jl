import Ogre
using FactCheck

facts("Heat capacity handling") do
    capfunc(T) = exp(T)

    Cₚ = Ogre.HeatCapacity(capfunc)

    @fact Cₚ(1) => roughly(e)
    @fact Cₚ(0) => 1
    @fact_throws Cₚ(1,2)

    vs1 = Ogre.ValueSet(1,2,3)
    vs2 = Ogre.ValueSet(1,2,3,4)

    @fact_throws Cₚ(vs1)
    @fact Cₚ(vs2) => exp(4)
end


