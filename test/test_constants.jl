if VERSION < v"0.5"
    using BaseTestNext
else
    using Base.Test
end

@testset "No constant tests" begin
    nothing
end
