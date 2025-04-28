using Test
using Jupo

@testset "base api" begin
    using DynamicalSystems

    @inbounds function lorenz_rule!(du, u, p, t)
        σ, r, b = p
        du[1] = σ * (u[2] - u[1])
        du[2] = - u[1] * u[3] + r * u[1] - u[2]
        du[3] = u[1] * u[2] - b * u[3]
        return nothing
    end

    function lorenz_rule_j!(du, u, p, t)
        σ, r, b = p
        du[1, :] = [-σ, σ, 0]
        du[2, :] = [(r - u[3]), -1, -u[1]]
        du[3, :] = [u[2], u[1], -b]
        return nothing
    end

    lorenz_test = System(
            "lrz63",
            lorenz_rule!,
            lorenz_rule_j!,
            3,
            [10., 28., 8/3],
            ["sig", "r", "b"]
        )

    @testset "system struct" begin
        @test lorenz_test isa System
        @test SystemDef(lorenz_test) isa SystemDef
    end

    @testset "TDS test" begin
        lrz_ode = generate_ode(lorenz_test)
        @test lrz_ode isa CoupledODEs

        lrz_ode, lrz_tds = generate_ds(lorenz_test)
        @test lrz_ode isa CoupledODEs
        @test lrz_tds isa TangentDynamicalSystem
    end
end

