using Test

# Cross-platform repo lint. Two kinds of checks:
#
#   * Path-literal lint: walk every .jl file under src/, test/, and
#     inputs/ and flag string literals that bake in platform-specific
#     path syntax — relative-parent segments (`../foo`, `..\foo`),
#     absolute Unix/Windows roots, Unix-only system roots (`/tmp`,
#     `/Users/...`), home-shorthand (`~/`), and `"$var/foo"` string
#     interpolation that builds a path with a hard-coded separator.
#     Such literals should be built with joinpath(...) instead so the
#     same code runs on macOS, Linux, and Windows.
#
#   * Filename lint: walk every tracked file in the repo and flag
#     basenames that Windows refuses to create — names containing any
#     of `< > : " \ | ? *`, or matching one of the reserved device
#     names (CON, PRN, AUX, NUL, COM1–COM9, LPT1–LPT9), with or without
#     an extension. A repo containing such a file can't be cloned on
#     Windows.

const _LINT_ROOT = normpath(joinpath(@__DIR__, ".."))
const _LINT_DIRS = ("src", "test", "inputs")

# Skip checked-in upstream reference material — it's read-only and may
# contain example paths that should not be rewritten.
const _LINT_SKIP_DIRS = (joinpath(_LINT_ROOT, "test", "simon-lab"),)

# Allow this file itself to mention the offending patterns in its
# docstrings and regex literals without self-flagging.
const _LINT_SKIP_FILES = (@__FILE__,)

"""
Strip the trailing line comment from `line`, respecting double-quoted
string literals (so a `#` inside `"foo#bar"` is not treated as a
comment). Returns the part of the line that is *not* a comment.
"""
function _strip_line_comment(line::AbstractString)
    io = IOBuffer()
    in_str = false
    chars = collect(line)
    i = 1
    n = length(chars)
    while i <= n
        c = chars[i]
        if in_str
            print(io, c)
            if c == '\\' && i < n
                print(io, chars[i + 1])
                i += 2
                continue
            elseif c == '"'
                in_str = false
            end
        else
            if c == '#'
                break
            elseif c == '"'
                in_str = true
                print(io, c)
            else
                print(io, c)
            end
        end
        i += 1
    end
    return String(take!(io))
end

# Unix-only filesystem roots. Each must be followed by `/` or end of
# string so that a math identifier like `/dev/null` triggers but a
# title like `"/dev"` does not... actually we want both. Use
# whole-prefix match on `"<root>/"` or `"<root>"`-then-end.
const _UNIX_ROOTS = ("tmp", "var", "etc", "usr", "opt", "home",
                     "private", "Users", "mnt", "dev", "root", "proc",
                     "sys", "bin", "sbin")

_unix_root_match(body) =
    any(occursin(Regex("^/" * r * raw"(/|$)"), body) for r in _UNIX_ROOTS)

"""
Return a vector of `(lineno, match, reason)` tuples for every
non-portable path literal in `path`.
"""
function lint_path_literals(path::AbstractString)
    issues = Tuple{Int, String, String}[]
    # Each rule is a `(predicate, reason)` pair applied to the body of
    # a double-quoted string literal (without the surrounding quotes).
    # Order matters: more specific rules first so the error message
    # points at the real cause.
    rules = [
        (body -> occursin("../",  body),
            "relative-parent segment '../' — use joinpath(\"..\", ...)"),
        (body -> occursin("..\\", body),
            "relative-parent segment '..\\\\' — use joinpath(\"..\", ...)"),
        (body -> startswith(body, "~/") || startswith(body, "~\\"),
            "tilde home shorthand — use `homedir()` (Julia does not expand `~`)"),
        (_unix_root_match,
            "hard-coded Unix-only root (e.g. /tmp, /Users) — build paths from a configurable base"),
        (body -> occursin(r"^/[A-Za-z]", body),
            "Unix-absolute path — build from joinpath / environment"),
        (body -> occursin(raw"^[A-Za-z]:[/\\]" |> Regex, body),
            "Windows-absolute path — build from joinpath / environment"),
        # Only flag interpolation-with-slash when the literal has no spaces;
        # status messages like "2/3 built" contain spaces and aren't paths.
        (body -> !occursin(' ', body) &&
                 (occursin(r"\$[A-Za-z_][A-Za-z0-9_]*/", body) ||
                  occursin(r"\$\([^)]*\)/", body)),
            "string-interpolated path with hard-coded '/' — use joinpath(\$var, ...)"),
    ]
    open(path) do io
        lineno = 0
        for line in eachline(io)
            lineno += 1
            stripped = _strip_line_comment(line)
            for m in eachmatch(r"\"((?:[^\"\\]|\\.)*)\"", stripped)
                body = m.captures[1]
                for (pred, reason) in rules
                    if pred(body)
                        push!(issues, (lineno, m.match, reason))
                        break
                    end
                end
            end
        end
    end
    return issues
end

function _lint_jl_targets()
    files = String[]
    # Top-level .jl files in the project root (BIFROST keeps source here,
    # not under src/). Non-recursive so we don't double-walk test/.
    for f in readdir(_LINT_ROOT)
        endswith(f, ".jl") || continue
        full = joinpath(_LINT_ROOT, f)
        isfile(full) || continue
        full in _LINT_SKIP_FILES && continue
        push!(files, full)
    end
    for d in _LINT_DIRS
        root = joinpath(_LINT_ROOT, d)
        isdir(root) || continue
        for (dir, _, fnames) in walkdir(root)
            any(startswith(dir, skip) for skip in _LINT_SKIP_DIRS) && continue
            for f in fnames
                endswith(f, ".jl") || continue
                full = joinpath(dir, f)
                full in _LINT_SKIP_FILES && continue
                push!(files, full)
            end
        end
    end
    return sort!(files)
end

# --- Filename lint -------------------------------------------------------

# Characters Windows rejects in any path component.
const _WINDOWS_FORBIDDEN_CHARS = ('<', '>', ':', '"', '\\', '|', '?', '*')

# Names Windows reserves regardless of extension. Match is
# case-insensitive against the component's stem.
const _WINDOWS_RESERVED_NAMES = Set([
    "CON", "PRN", "AUX", "NUL",
    "COM1", "COM2", "COM3", "COM4", "COM5", "COM6", "COM7", "COM8", "COM9",
    "LPT1", "LPT2", "LPT3", "LPT4", "LPT5", "LPT6", "LPT7", "LPT8", "LPT9",
])

"""
Return reasons (possibly multiple) why `relpath` would be invalid on
Windows. `relpath` is the repo-relative path with `/` separators.
"""
function lint_filename(relpath::AbstractString)
    reasons = String[]
    for component in split(relpath, '/'; keepempty = false)
        for c in _WINDOWS_FORBIDDEN_CHARS
            if c in component
                push!(reasons,
                      "component '$component' contains '$c' (Windows-forbidden in filenames)")
                break
            end
        end
        stem = first(split(component, '.'; limit = 2))
        if uppercase(stem) in _WINDOWS_RESERVED_NAMES
            push!(reasons,
                  "component '$component' is a reserved Windows device name ($(uppercase(stem)))")
        end
    end
    return reasons
end

function _tracked_files()
    # Prefer `git ls-files`: it's the authoritative shipped-file list,
    # and it already skips .git/, outputs/, and other ignored paths.
    out = try
        read(Cmd(`git ls-files`; dir = _LINT_ROOT), String)
    catch
        ""
    end
    if !isempty(out)
        return [strip(l) for l in split(out, '\n') if !isempty(strip(l))]
    end
    # Fallback: walk the tree skipping a hard-coded ignore list.
    skip_dirs = ("/.git", "/outputs", "/scratch", "/.DS_Store")
    files = String[]
    for (dir, _, fnames) in walkdir(_LINT_ROOT)
        rel_dir = "/" * relpath(dir, _LINT_ROOT)
        any(occursin(s, rel_dir) for s in skip_dirs) && continue
        for f in fnames
            push!(files, relpath(joinpath(dir, f), _LINT_ROOT))
        end
    end
    return files
end

# --- Tests ---------------------------------------------------------------

@testset "cross-platform lint" begin
    @testset "path literals" begin
        files = _lint_jl_targets()
        @test !isempty(files)
        offenders = Tuple{String, Int, String, String}[]
        for f in files
            for (ln, lit, reason) in lint_path_literals(f)
                push!(offenders, (relpath(f, _LINT_ROOT), ln, lit, reason))
            end
        end
        if !isempty(offenders)
            msg = IOBuffer()
            println(msg, "Non-portable path literals detected:")
            for (rel, ln, lit, reason) in offenders
                println(msg, "  $rel:$ln  $lit  — $reason")
            end
            @error String(take!(msg))
        end
        @test isempty(offenders)
    end

    @testset "filenames" begin
        files = _tracked_files()
        @test !isempty(files)
        offenders = Tuple{String, String}[]
        for f in files
            for reason in lint_filename(f)
                push!(offenders, (f, reason))
            end
        end
        if !isempty(offenders)
            msg = IOBuffer()
            println(msg, "Windows-hostile filenames detected:")
            for (f, reason) in offenders
                println(msg, "  $f — $reason")
            end
            @error String(take!(msg))
        end
        @test isempty(offenders)
    end

    @testset "linter self-check" begin
        # Path-literal rules: each forbidden pattern is caught, and
        # legitimate joinpath / comment forms are not.
        mktempdir() do dir
            sample = joinpath(dir, "sample.jl")
            open(sample, "w") do io
                println(io, "a = \"../outputs/foo.csv\"")           # 1: ../
                println(io, "b = \"..\\\\outputs\\\\foo.csv\"")     # 2: ..\
                println(io, "c = \"~/data/foo.csv\"")               # 3: ~/
                println(io, "d = \"/tmp/foo.csv\"")                 # 4: unix root /tmp
                println(io, "e = \"/Users/jane/foo.csv\"")          # 5: unix root /Users
                println(io, "f = \"/opt/local/share\"")             # 6: unix root /opt
                println(io, "g = \"C:\\\\Users\\\\jane\\\\f.csv\"") # 7: windows-absolute
                println(io, "h = \"\$dir/file.csv\"")               # 8: interp + /
                println(io, "i = \"\$(cfg.outdir)/file.csv\"")      # 9: interp(...) + /
                println(io, "ok1 = joinpath(\"..\", \"outputs\", \"f.csv\")")
                println(io, "ok2 = \"\$(cfg.name)-fir.csv\"")
                println(io, "# comment with ../foo and /tmp/bar is ignored")
            end
            hits = lint_path_literals(sample)
            lines = Set(h[1] for h in hits)
            @test issetequal(lines, Set(1:9))
        end

        # Filename rules: forbidden chars and reserved names are caught.
        @test isempty(lint_filename("src/normal_file.jl"))
        @test isempty(lint_filename("inputs/20260513a_M1.csv"))
        @test !isempty(lint_filename("src/foo:bar.jl"))
        @test !isempty(lint_filename("src/what?.jl"))
        @test !isempty(lint_filename("dir/CON.txt"))
        @test !isempty(lint_filename("dir/nul"))
        @test !isempty(lint_filename("LPT3.csv"))
        # Reserved name buried inside a longer word is OK.
        @test isempty(lint_filename("inputs/console.log"))
    end
end
