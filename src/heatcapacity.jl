# Evaluating heat capacities

using WaterData
# needed so that WaterData is able to load files
import JLD, Dierckx, VoronoiDelaunay, GeometricalPredicates

# add an override to call heat capacities with ValueSets
for HC in (WaterData.PTFuncHeatCapacity, WaterData.GridHeatCapacity,
           WaterData.TFuncHeatCapacity, WaterData.ConstantHeatCapacity)
    Base.call(cp::HC, pv::PhysicalValues) = cp(pressure(pv), temperature(pv))
end
