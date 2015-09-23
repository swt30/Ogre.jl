using FactCheck
import Ogre


module test_heatcapacity_resources
    import Ogre

    func = exp
    exponential = Ogre.TFuncHeatCapacity(func)
    func2d(x, y) = exp(x) + 2*exp(y)
    exponential2d = Ogre.PTFuncHeatCapacity(func2d)
    constant = Ogre.ConstantHeatCapacity(1)
end


facts("Heat capacity handling") do
    res = test_heatcapacity_resources

    context("Values are as expected") do
        cₚ_exp = res.exponential
        cₚ_exp2 = res.exponential2d
        cₚ_const = res.constant

        context("called directly") do
            # constant values are always equal
            @fact cₚ_const(0) --> cₚ_const(10)
            @fact cₚ_const(0, 10) --> 1

            # calling a single-input heat capacity
            @fact cₚ_exp(1) --> e^1
            @fact cₚ_exp(0) --> e^0
            @fact cₚ_exp(1, 2) --> e^2

            # calling a temp- and pressure-dependent heat capacity
            @fact cₚ_exp2(0, 1) --> e^0 + 2e^1
            @fact cₚ_exp2(1, 1) --> e^1 + 2e^1
            @fact cₚ_exp2(1, 3) --> e^1 + 2e^3
            @fact cₚ_exp2(2, 1) --> e^2 + 2e^1

            # needs both temperature and pressure
            @fact_throws MethodError cₚ_exp2(1)
        end

        context("called through a ValueSet") do
            vs_notemp = Ogre.ValueSet(1,2,3)
            vs_full = Ogre.ValueSet(1,2,3,4)
            # accessing temperature and pressure
            @fact_throws MethodError cₚ_exp(vs_notemp)
            @fact Ogre.temperature(vs_full) --> 4
            @fact Ogre.pressure(vs_full) --> 3

            # uses just temperature
            @fact cₚ_exp(vs_full) --> e^4

            # uses temperature and pressure
            @fact cₚ_exp2(vs_full) --> e^3 + 2e^4
        end
    end
end
