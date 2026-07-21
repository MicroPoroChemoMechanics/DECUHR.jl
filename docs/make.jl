using Documenter
using Decuhr

DocMeta.setdocmeta!(
    Decuhr,
    :DocTestSetup,
    :(using Decuhr);
    recursive = true,
)

makedocs(
    clean    = false,
    modules  = [Decuhr],
    remotes  = nothing,
    authors  = "Jean-François Barthélémy",
    sitename = "Decuhr.jl",
    format   = Documenter.HTML(;
        canonical     = "https://MicroPoroChemoMechanics.github.io/Decuhr.jl",
        repolink      = "https://github.com/MicroPoroChemoMechanics/Decuhr.jl",
        edit_link     = "main",
        assets        = ["assets/custom.css"],
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
    repo         = "github.com/MicroPoroChemoMechanics/Decuhr.jl.git",
    devbranch    = "main",
    push_preview = false,
)
