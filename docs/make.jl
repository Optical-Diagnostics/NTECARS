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
        "References" => [
            "Spectra"  => "reference/spectra.md",
            "CARSSimulator setup"  => "reference/CARSSimulator_construction.md",
            "CO2" => "reference/CO2.md",
            "N2" => "reference/N2.md",
            "Fitting"  => "reference/fitting.md"
        ]#"reference.md"
    ],
    #clean=true,
    checkdocs=:exports,
    #remotes = nothing,
)

deploydocs(
    repo = "github.com/Optical-Diagnostics/NTECARS.git",
)