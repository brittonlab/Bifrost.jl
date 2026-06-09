using Test

const _FILE_TESTSETS = Test.AbstractTestSet[]

function _run_test_file(path)
    label = splitext(basename(path))[1]
    ts = @testset "$label" begin
        include(path)
    end
    push!(_FILE_TESTSETS, ts)
end

_run_test_file(joinpath(@__DIR__, "test_path_geometry_connector.jl"))
_run_test_file(joinpath(@__DIR__, "test_path_geometry.jl"))
_run_test_file(joinpath(@__DIR__, "test_path_geometry_perturb.jl"))
_run_test_file(joinpath(@__DIR__, "test_fiber_thermal.jl"))
_run_test_file(joinpath(@__DIR__, "test_fiber_tension.jl"))
_run_test_file(joinpath(@__DIR__, "test_fiber_path_pass3.jl"))
_run_test_file(joinpath(@__DIR__, "test_material_properties.jl"))
_run_test_file(joinpath(@__DIR__, "test_mcm_compatability.jl"))
_run_test_file(joinpath(@__DIR__, "test_paddle_transfer.jl"))
_run_test_file(joinpath(@__DIR__, "test_dgd.jl"))
_run_test_file(joinpath(@__DIR__, "test_fiber_cross_section.jl"))
_run_test_file(joinpath(@__DIR__, "test_path_integral.jl"))
_run_test_file(joinpath(@__DIR__, "test_solver_convergence.jl"))
_run_test_file(joinpath(@__DIR__, "cross-platform_tests.jl"))

let passes = 0, fails = 0, errors = 0
    for ts in _FILE_TESTSETS
        c = Test.get_test_counts(ts)
        passes += c.passes + c.cumulative_passes
        fails  += c.fails  + c.cumulative_fails
        errors += c.errors + c.cumulative_errors
    end
    println("BIFROST tests: $passes passed, $fails failed, $errors errored")
end
