using Documenter, DocumenterCitations, Literate

using Walrus
#=
using CairoMakie
CairoMakie.activate!(type = "svg")
=#
bib_filepath = joinpath(dirname(@__FILE__), "bibliography.bib")
bib = CitationBibliography(bib_filepath, style=:authoryear)

# Examples

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/generated")

examples = [
    #"Single plant" => "single_plant",
    #"Forest" => "forest"
]

example_scripts = [ filename * ".jl" for (title, filename) in examples ]

for example in example_scripts
    example_filepath = joinpath(EXAMPLES_DIR, example)

    withenv("JULIA_DEBUG" => "Literate") do
        Literate.markdown(example_filepath, OUTPUT_DIR; 
                          flavor = Literate.DocumenterFlavor(),
                          repo_root_url = "https://jagows.com/Walrus.jl",
                          execute = true)
    end
end

example_pages = [ title => "generated/$(filename).md" for (title, filename) in examples ]

numerical_pages = [
    "Coming soon" => "coming-soon.md"
]

appendix_pages = [
    "Library" => "appendix/library.md",
    "Function index" => "appendix/function_index.md",
]

pages = [
    "Home" => "index.md",
    "Quick start" => "coming-soon.md",
    "Examples" => example_pages,
    "Numerical implementation" => numerical_pages,
    "References" => "references.md",
    "Appendix" => appendix_pages
]

#####
##### Build and deploy docs
#####

format = Documenter.HTML(
    collapselevel = 1,
    prettyurls = get(ENV, "CI", nothing) == "true",
    canonical = "https://jagosw.com/Walrus/stable/",
    mathengine = MathJax3(),
    assets = String["assets/citations.css"]
)

makedocs(sitename = "Walrus.jl",
         authors = "Jago Strong-Wright",
         format = format,
         pages = pages,
         modules = [Walrus],
         plugins = [bib],
         doctest = true,
         clean = true,
         checkdocs = :exports)

@info "Clean up temporary .jld2/.nc files created by doctests..."

"""
    recursive_find(directory, pattern)

Return list of filepaths within `directory` that contains the `pattern::Regex`.
"""
recursive_find(directory, pattern) =
    mapreduce(vcat, walkdir(directory)) do (root, dirs, files)
        joinpath.(root, filter(contains(pattern), files))
    end

files = []
for pattern in [r"\.jld2", r"\.nc"]
    global files = vcat(files, recursive_find(@__DIR__, pattern))
end

for file in files
    rm(file)
end

deploydocs(
    repo = "github.com/jagoosw/GiantKelpDynamics",
    versions = ["stable" => "v^", "dev" => "dev", "v#.#.#"],
    forcepush = true,
    push_preview = true,
    devbranch = "main"
)
