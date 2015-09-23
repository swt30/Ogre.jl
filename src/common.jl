# Types which indicate whether the system includes temperature dependence


"Define the complexity of a physical system"
abstract ModelComplexity
"This system excludes temperature details"
immutable NoTemp <: ModelComplexity; end
"This system explicitly includes temperature details"
immutable WithTemp <: ModelComplexity; end
"This system explicitly includes both temperature and pressure details"
immutable WithTempPressure <: ModelComplexity; end
