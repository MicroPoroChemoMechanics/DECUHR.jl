using Documenter
using Documenter.Remotes
using Pkg

# Activate the package so DECUHR is loadable from the docs environment
Pkg.develop(PackageSpec(path = joinpath(@__DIR__, "..")))

using DECUHR

DocMeta.setdocmeta!(
    DECUHR,
    :DocTestSetup,
    :(using DECUHR);
    recursive = true,
)

makedocs(
    clean    = false,
    modules  = [DECUHR],
    authors  = "Jean-François Barthélémy",
    sitename = "DECUHR.jl",
    remotes  = Dict(
        joinpath(@__DIR__, "..") => (Remotes.GitHub("MicMacTools", "DECUHR.jl"), "main"),
    ),
    format   = Documenter.HTML(;
        canonical        = "https://MicMacTools.github.io/DECUHR.jl",
        edit_link        = "main",
        assets           = ["assets/favicon.ico", "assets/custom.css"],
        prettyurls       = (get(ENV, "CI", nothing) == "true"),
        collapselevel    = 1,
        ansicolor        = true,
    ),
    pages = [
        "Home"          => "index.md",
        "Algorithm"     => "algorithm.md",
        "Examples"      => "examples.md",
        "API Reference" => "api.md",
        "License"       => "license.md",
    ],
    checkdocs = :exports,
    warnonly  = true,
)

deploydocs(;
    repo      = "github.com/MicMacTools/DECUHR.jl.git",
    devbranch = "main",
)
