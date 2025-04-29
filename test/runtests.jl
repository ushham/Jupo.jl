using Test
using Jupo


@testset "code quality" begin
    include("quality_check.jl")
end

@testset "Base" begin
    include("base.jl")
end

@testset "Newton Method" begin
    include("newton_method.jl")
end

### SYSTEM CHECKS ###
@testset "lorenz system" begin
    include("system_tests/lorenztest.jl")
end

@testset "Land Atmosphere system" begin
    include("system_tests/maolamtest.jl")
end
