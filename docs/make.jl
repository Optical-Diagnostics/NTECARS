using Documenter
using LiveServer
using NTECARS  # load your package
#DocMeta.setdocmeta!(NTECARS, :DocTestSetup, :(using NTECARS); recursive=true)
makedocs(
    modules=[NTECARS],
    sitename="NTECARS Documentation",
    pages=[
        "Home" => "index.md",
        "Tutorial" => [
            "Setup for Spectra Calculation" => "Setup_for_spectra.md",
            "Fitting" => "fitting.md"
        ],
        "References" => "reference.md"
    ],
    #clean=true,
    checkdocs=:exports,
    #remotes = nothing,
)

deploydocs(
    repo = "github.com/ChristianABusch/NTECARS.git",
)