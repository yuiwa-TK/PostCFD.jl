using PostCFD
using Documenter

DocMeta.setdocmeta!(PostCFD, :DocTestSetup, :(using PostCFD); recursive=true)

makedocs(;
    modules=[PostCFD],
    authors="Yuta Iwatani <yuta.iwatani.cfd@gmail.com> and contributors",
    repo="https://github.com/yuiwa-TK/PostCFD.jl/blob/{commit}{path}#{line}",
    sitename="PostCFD.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://yuiwa-TK.github.io/PostCFD.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API Reference" => Any[
            "FileReader" => "filereader.md",
            "Filewriter" => "filewriter.md",
            "MathLib" => "mathlib.md",
        ]
    ],
)

deploydocs(;
    repo="github.com/yuiwa-TK/PostCFD.jl",
    devbranch="main",
)