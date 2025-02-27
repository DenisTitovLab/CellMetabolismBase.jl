using CellMetabolismBase
using Documenter

DocMeta.setdocmeta!(CellMetabolismBase, :DocTestSetup, :(using CellMetabolismBase); recursive=true)

makedocs(;
    modules=[CellMetabolismBase],
    authors="Denis Titov <titov@berkeley.edu>  and contributors",
    sitename="CellMetabolismBase.jl",
    format=Documenter.HTML(;
        canonical="https://Denis-Titov.github.io/CellMetabolismBase.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Denis-Titov/CellMetabolismBase.jl",
    devbranch="main",
)
