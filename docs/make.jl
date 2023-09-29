using DecoratedParticles
using Documenter

DocMeta.setdocmeta!(DecoratedParticles, :DocTestSetup, :(using DecoratedParticles); recursive=true)

makedocs(;
    modules=[DecoratedParticles],
    authors="Christoph Ortner <christohortner@gmail.com> and contributors",
    repo="https://github.com/ACEsuit/DecoratedParticles.jl/blob/{commit}{path}#{line}",
    sitename="DecoratedParticles.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ACEsuit.github.io/DecoratedParticles.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ACEsuit/DecoratedParticles.jl",
    devbranch="main",
)
