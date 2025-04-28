using Test
using Jupo

@testset "quality_check" begin
    using Aqua, Documenter
    ignore_deps = [:LinearAlgebra, :Test, :Plots]
    Aqua.test_all(Jupo;
        deps_compat=(
                ignore=ignore_deps,
                check_extras=(ignore=ignore_deps,),
                check_weakdeps=(ignore=ignore_deps,),
            )
    )
end

# Documentation not complete

# @testset "Docs test" begin
#     Documenter.doctest(Jupo)
# end