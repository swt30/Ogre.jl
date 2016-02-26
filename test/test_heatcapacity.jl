using FactCheck
import Ogre


module test_heatcapacity_resources
import Ogre
import WaterData

func = exp
exponential = WaterData.TFuncHeatCapacity(func)
func2d(x, y) = exp(x) + 2*exp(y)
exponential2d = WaterData.PTFuncHeatCapacity(func2d)
constant = WaterData.ConstantHeatCapacity(1)
end


facts("Heat capacity handling") do
    res = test_heatcapacity_resources
    context("called through a ValueSet") do
        vs_notemp = Ogre.ValueSet(1,2,3)
        vs_full = Ogre.ValueSet(1,2,3,4)
        # accessing temperature and pressure
        @fact_throws MethodError res.exponential(vs_notemp)
        @fact Ogre.temperature(vs_full) --> 4
        @fact Ogre.pressure(vs_full) --> 3

        # uses just temperature
        @fact res.exponential(vs_full) --> e^4

        # uses temperature and pressure
        @fact res.exponential2d(vs_full) --> e^3 + 2e^4
    end
end
