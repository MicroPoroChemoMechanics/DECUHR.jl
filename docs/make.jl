using Documenter
using DECUHR

# Generates the Examples page (plus a companion notebook and cleaned script)
# from `examples/basic_usage.jl` before `makedocs` runs, so the generated
# markdown exists when `pages` below references it.
#
# Before paying for a full build, the cheap guard that catches the one mistake
# able to break every page at once — a `Pkg` call that lost its `#jl` marker:
#
#     grep -n 'Pkg\.' examples/basic_usage.jl | grep -v '#jl'   # prints nothing
#
# `docs/literate.jl` re-checks it at build time and errors out early.
include("literate.jl")

DocMeta.setdocmeta!(
    DECUHR,
    :DocTestSetup,
    :(using DECUHR);
    recursive = true,
)

makedocs(
    clean    = false,
    modules  = [DECUHR],
    remotes  = nothing,
    authors  = "Jean-François Barthélémy",
    sitename = "DECUHR.jl",
    format   = Documenter.HTML(;
        canonical     = "https://MicroPoroChemoMechanics.github.io/DECUHR.jl",
        repolink      = "https://github.com/MicroPoroChemoMechanics/DECUHR.jl",
        edit_link     = "main",
        assets        = ["assets/custom.css"],
        prettyurls    = (get(ENV, "CI", nothing) == "true"),
        collapselevel = 1,
        ansicolor     = true,
    ),
    pages = [
        "Home"          => "index.md",
        "Algorithm"     => "algorithm.md",
        "Examples"      => "generated/examples.md",
        "API Reference" => "api.md",
        "License"       => "license.md",
    ],
    checkdocs = :exports,
    warnonly  = true,
)

deploydocs(;
    repo         = "github.com/MicroPoroChemoMechanics/DECUHR.jl.git",
    devbranch    = "main",
    push_preview = false,
)
