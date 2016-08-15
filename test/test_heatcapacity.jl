using Base.Test
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


@testset "Heat capacity handling" begin
    res = test_heatcapacity_resources
    @testset "called through a ValueSet" begin
        vs_notemp = Ogre.ValueSet(1,2,3)
        vs_full = Ogre.ValueSet(1,2,3,4)
        # accessing temperature and pressure
        @test_throws MethodError res.exponential(vs_notemp)
        @test Ogre.temperature(vs_full) == 4
        @test Ogre.pressure(vs_full) == 3

        # uses just temperature
        @test res.exponential(vs_full) == e^4

        # uses temperature and pressure
        @test res.exponential2d(vs_full) == e^3 + 2e^4
    end
end
