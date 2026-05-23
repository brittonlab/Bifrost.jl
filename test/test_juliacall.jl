using Test

# Drive the Python → Julia (juliacall) probes from the standard Julia test
# infrastructure. Each probe is a standalone Python script that boots Julia
# via `bifrost.start()` and exercises a specific juliacall failure mode that
# the native Julia test suite cannot reach (Julia version selection,
# multithreading, signal handling, project environments, MCM under threads,
# numerical parity between juliacall and native Julia).
#
# Probes are launched as subprocesses with `uv run python <script>`, so each
# gets a clean Julia runtime — that's the whole point of these tests.

const _JULIACALL_REPO = normpath(joinpath(@__DIR__, ".."))

const _JULIACALL_PROBES = (
    joinpath(_JULIACALL_REPO, "docs", "juliacall-demo.py"),
    joinpath(_JULIACALL_REPO, "test", "juliacall", "juliacall-mcm.py"),
)

function _run_python_probe(script::AbstractString)
    cmd = Cmd(`uv run python $script`; dir = _JULIACALL_REPO)
    return success(pipeline(cmd; stdout = stdout, stderr = stderr))
end

@testset "juliacall" begin
    # If `uv` isn't on PATH we can't drive these probes — treat that as a
    # broken environment rather than a silent skip, but make the failure
    # message actionable.
    uv = Sys.which("uv")
    @test uv !== nothing

    if uv !== nothing
        for script in _JULIACALL_PROBES
            name = relpath(script, _JULIACALL_REPO)
            @testset "$name" begin
                @test isfile(script)
                @test _run_python_probe(script)
            end
        end
    end
end
