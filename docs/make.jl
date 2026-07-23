using Documenter, SpectralKurtosisEstimators

makedocs(
    sitename = "SpectralKurtosisEstimators.jl",
    authors = "David MacMahon",
    modules = [SpectralKurtosisEstimators],
    format = Documenter.HTML(
        canonical = "https://david-macmahon.github.io/SpectralKurtosisEstimators.jl/stable/",
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "SKEstimator" => "ske.md",
        "Pearson Distributions" => "pearson_distributions.md",
    ],
    doctest = true,
    checkdocs = :exports,
    remotes = nothing,
)

deploydocs(
    repo = "github.com/david-macmahon/SpectralKurtosisEstimators.jl.git",
    devbranch = "main",
)
