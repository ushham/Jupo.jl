using Test
using Jupo


# multiple dispatch part failing
# @testset "code quality" begin
#     include("quality_check.jl")
# end


@testset "Base" begin
    include("base.jl")
end

# @testset "Newton Method" begin
#     include("newton_method.jl")
# end

### SYSTEM CHECKS ###
# @testset "lorenz system" begin
#     include("system_tests/lorenztest.jl")
# end
