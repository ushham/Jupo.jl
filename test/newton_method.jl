using Test
using Jupo

@testset "Newton Method" begin
    using DynamicalSystems
    using Jupo:
        construct_matrix!

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

    @testset "Matrix construction" begin
        lrz_ode, lrz_tds = generate_ds(lorenz_test)

        x_i = rand(3)
        vec_1 = rand(3)
        lhs_mat = zeros(3, 3)
        rhs_vec = rand(3)
        time = 1.
        vec_2 = rand(3)
        vec_3 = rand(3)

        @test construct_matrix!(
            tds=lrz_tds, 
            x_i = x_i,
            t=1.,
            ndim=3,
            s_x=vec_1,
            lhs_mat=lhs_mat,
            rhs_vec=rhs_vec,
            f_x=vec_2,
            f_s_x=vec_3,
            f=lorenz_test.f
        ) isa Nothing
    end

    # Tests on tensor section to be added

end

