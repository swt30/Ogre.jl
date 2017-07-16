# Types which indicate whether the system includes temperature dependence


"Define the complexity of a physical system"
abstract type ModelComplexity end
"This system excludes temperature details"
struct NoTemp <: ModelComplexity end
"This system explicitly includes temperature details"
struct WithTemp <: ModelComplexity end
"This system explicitly includes both temperature and pressure details"
struct WithTempPressure <: ModelComplexity end
