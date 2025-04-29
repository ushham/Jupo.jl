using Test
using Jupo

@testset "Lorenz63 System" begin
    
    @testset "sys" begin
        @test lorenz isa System
    end

    @testset "upo NM" begin
        x_ic = [14.1, 9.5, 38.64, 1.55]
        x = find_upo_nm(lorenz, x_ic)

        @test x isa UPO_sol
        @test x.period ≈ 1.5586522107124108
    end
end