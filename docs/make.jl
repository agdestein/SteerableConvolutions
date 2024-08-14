# Look for environment variable triggering local development modifications
localdev = haskey(ENV, "LOCALDEV")

# Get access to example dependencies
push!(LOAD_PATH, joinpath(@__DIR__, "..", "examples"))

using SteerableConvolutions
using Literate
using Documenter
using DocumenterCitations
using DocumenterVitepress

DocMeta.setdocmeta!(
    SteerableConvolutions,
    :DocTestSetup,
    :(using SteerableConvolutions);
    recursive = true,
)

bib = CitationBibliography(joinpath(@__DIR__, "references.bib"))

# Generate examples
examples = [
    (true, "examples/introduction"),
    (true, "examples/equivariance"),
    # (true, "examples/scratch"),
]

# Convert scripts to executable markdown files
output = "examples/generated"
outputdir = joinpath(@__DIR__, "src", output)
## rm(outputdir; recursive = true)
for (run, name) in examples
    inputfile = joinpath(@__DIR__, "..", name * ".jl")
    if run
        # With code execution blocks
        Literate.markdown(inputfile, outputdir)
    else
        # Turn off code execution.
        # Note: Literate has a `documenter = false` option, but this would also remove
        # the "Edit on GitHub" button at the top, therefore we disable the `@example`-blocks
        # manually
        Literate.markdown(
            inputfile,
            outputdir;
            preprocess = content ->
                "# *Note: Output is not generated for this example (to save resources on GitHub).*\n\n" *
                content,
            postprocess = content -> replace(content, r"@example.*" => "julia"),
        )
    end
end

vitepress_kwargs = localdev ? (;
    # md_output_path = @__DIR__,
    build_vitepress = false
) : (;)

makedocs(;
    # draft = true,
    # clean = !localdev,
    modules = [SteerableConvolutions],
    plugins = [bib],
    authors = "Syver Døving Agdestein and contributors",
    repo = Remotes.GitHub("agdestein", "SteerableConvolutions.jl"),
    sitename = "SteerableConvolutions.jl",
    format = DocumenterVitepress.MarkdownVitepress(;
        repo = "github.com/agdestein/SteerableConvolutions.jl",
        devurl = "dev",
        vitepress_kwargs...,
    ),
    pagesonly = true,
)

# Only deploy docs on CI
get(ENV, "CI", "false") == "true" && deploydocs(;
    repo = "github.com/agdestein/SteerableConvolutions.jl",
    target = "build",
    devbranch = "main",
    push_preview = true,
)
