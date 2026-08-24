# docs/literate.jl
#
# Runs Literate.jl over `examples/` before `makedocs`. The example file produces
# three artifacts from a single source:
#
#   - a Documenter-ready markdown page  -> docs/src/generated/
#     (executed by Documenter itself through `@example`; Literate's
#     `markdown(...; documenter = true)` leaves `execute = false`)
#   - a pre-run Jupyter notebook        -> docs/generated_notebooks/
#   - a cleaned standalone .jl script   -> docs/generated_scripts/
#
# `docs/make.jl` includes this file *before* `makedocs`, so the generated
# markdown exists by the time the `pages` list references it.
#
# The example stays runnable on its own (`julia --project=. examples/basic_usage.jl`);
# the Literate markers are the only thing that distinguishes it from an ordinary
# script. Before this file existed, `examples/basic_usage.jl` and
# `docs/src/examples.md` were two hand-maintained copies of the same integrals,
# and the page had quietly grown three sections the script never gained.
#
# ┌─────────────────────────────────────────────────────────────────────────┐
# │ The `#jl` contract, which is not optional                               │
# │                                                                         │
# │ Every `import Pkg` / `Pkg.activate` / `Pkg.instantiate` line MUST end    │
# │ with `#jl`. That marker keeps the line in the standalone script and in   │
# │ the generated clean script, and strips it from the markdown and the      │
# │ notebook — which is what we want, because those already run inside the   │
# │ `docs` environment.                                                     │
# │                                                                         │
# │ Forgetting it does not break only that page. `Pkg.activate` mutates      │
# │ global process state, so Documenter switches project mid-build and every │
# │ later `@example` block, on every page, fails with                        │
# │ `Package DECUHR not found in current path`. `check_pkg_markers()` below  │
# │ turns that rule into an assertion, checked on every build.               │
# └─────────────────────────────────────────────────────────────────────────┘

using Literate

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const GENERATED_MD_DIR = joinpath(@__DIR__, "src", "generated")
const NOTEBOOK_DIR = joinpath(@__DIR__, "generated_notebooks")
const CLEAN_SCRIPT_DIR = joinpath(@__DIR__, "generated_scripts")

# Source file => page name.
const PUBLISHED_EXAMPLES = ["basic_usage.jl" => "examples"]

"""
    check_pkg_markers()

Fail loudly if a published example carries a bare `Pkg.activate` — the one
mistake that breaks the whole build rather than a single page. Cheap enough to
run on every build.
"""
function check_pkg_markers()
    offenders = String[]
    for (script, _) in PUBLISHED_EXAMPLES
        path = joinpath(EXAMPLES_DIR, script)
        isfile(path) || continue
        for line in eachline(path)
            occursin(r"^\s*(import\s+Pkg|Pkg\.(activate|instantiate))", line) &&
                !occursin("#jl", line) &&
                push!(offenders, "$script: $(strip(line))")
        end
    end
    isempty(offenders) || error(
        "Published examples carry a bare `Pkg` call (missing the `#jl` marker).\n" *
            "This would switch the active project mid-build and break every @example\n" *
            "block in the documentation, not just these pages:\n  " *
            join(offenders, "\n  ")
    )
    return nothing
end

function build_example_pages()
    check_pkg_markers()
    mkpath(GENERATED_MD_DIR)
    mkpath(NOTEBOOK_DIR)
    mkpath(CLEAN_SCRIPT_DIR)
    for (script, page) in PUBLISHED_EXAMPLES
        src = joinpath(EXAMPLES_DIR, script)
        isfile(src) || error("docs/literate.jl: published example not found: $src")
        Literate.markdown(src, GENERATED_MD_DIR; documenter = true, name = page)
        Literate.notebook(src, NOTEBOOK_DIR; name = page)
        Literate.script(src, CLEAN_SCRIPT_DIR; name = page)
    end
    return nothing
end

build_example_pages()
