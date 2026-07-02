using Documenter
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
    remotes  = nothing,
    authors  = "Jean-François Barthélémy",
    sitename = "DECUHR.jl",
    format   = Documenter.HTML(;
        canonical     = "https://MicroPoroChemoMechanics.github.io/DECUHR.jl",
        repolink      = "https://github.com/MicroPoroChemoMechanics/DECUHR.jl",
        edit_link     = "main",
        assets        = ["assets/favicon.ico", "assets/custom.css"],
        prettyurls    = (get(ENV, "CI", nothing) == "true"),
        collapselevel = 1,
        ansicolor     = true,
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
    repo         = "github.com/MicroPoroChemoMechanics/DECUHR.jl.git",
    devbranch    = "main",
    push_preview = false,
)
