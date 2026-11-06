using Documenter
using DECUHR
using DocumenterCitations
using DocumenterVitepress

# Generates the Examples page (plus a companion notebook and cleaned script)
# from `examples/basic_usage.jl` before `makedocs` runs, so the generated
# markdown exists when `pages` below references it.
#
# Before paying for a full build, the cheap guard that catches the one mistake
# able to break every page at once — a `Pkg` call that lost its `#jl` marker:
#
#     grep -n 'Pkg\.' examples/basic_usage.jl | grep -v '#jl'   # prints nothing
#
# `docs/literate.jl` re-checks it at build time and errors out early.
include("literate.jl")

DocMeta.setdocmeta!(
    DECUHR,
    :DocTestSetup,
    :(using DECUHR);
    recursive = true,
)

# Author-year citations, as in TensND, MeanFieldHomogenization and ChemistryLab:
# a `[3]` tells the reader nothing without a trip to the References page, and
# hovering a citation shows the full entry (`theme/index.ts`).
bib = CitationBibliography(
    joinpath(@__DIR__, "src", "references.bib");
    style = :authoryear,
)

# ── The bibliography is formatted NOW ────────────────────────────────────────
#
# DocumenterCitations formats the entries only at its last stage, once every
# page has been expanded and every example run. A malformed entry would then
# kill the build after all that work. Calling exactly what the late stage calls
# moves that failure to the first seconds. A LaTeX escape the TeX parser does not
# implement is the usual cause: write the character bare inside braces instead.
let failures = String[]
    for (key, entry) in bib.entries
        try
            DocumenterCitations.format_bibliography_reference(:authoryear, entry)
        catch err
            push!(failures, "  $key : " * sprint(showerror, err))
        end
    end
    isempty(failures) || error(
        "docs/src/references.bib has $(length(failures)) entry/entries " *
            "DocumenterCitations cannot format. This would otherwise kill the " *
            "build at its LAST stage:\n" *
            join(failures, "\n")
    )
end

# ── Stopgap: CitationSiteNode in the Markdown writer ─────────────────────────
#
# GUARDED, and the guard is what makes one file serve both versions.
#
# `CitationSiteNode` exists only in DocumenterCitations 1.5, while
# DocumenterVitepress declares a weak dependency `DocumenterCitations = "1 - 1.4"`
# — an upstream cap, so the docs environment resolves to 1.4 and the name is not
# there. Referring to it unconditionally is an `UndefVarError` raised before the
# first page is built. On 1.4 the citations need no help at all: the node is what
# 1.5 introduced.
#
# When DocumenterVitepress raises its cap, widening the bound in
# `docs/Project.toml` is the only edit needed — this method starts applying again
# on its own.
if isdefined(DocumenterCitations, :CitationSiteNode)
# DocumenterCitations 1.5 wraps every expanded citation in a `CitationSiteNode`,
# whose only purpose is to give the citation an HTML anchor so the bibliography
# can link back to it. Its own docstring calls it "transparent in any output
# format other than HTML", and both the LaTeX writer and MDFlatten implement it
# as "render my children".
#
# DocumenterVitepress 0.3.5 ships a DocumenterCitations extension, but it covers
# only `BibliographyNode`. With no method for `CitationSiteNode`, the writer
# falls through to its generic branch, which prints `Markdown.plain(element)` —
# so every citation would come out as the literal text
# `DocumenterCitations.CitationSiteNode("<key>-cite-1")`.
#
# The same one-line treatment as the other non-HTML writers. Remove this once
# DocumenterVitepress covers the node upstream.
    function DocumenterVitepress.render(
            io::IO,
            mime::MIME"text/plain",
            node::Documenter.MarkdownAST.Node,
            ::DocumenterCitations.CitationSiteNode,
            page,
            doc;
            kwargs...,
        )
        return DocumenterVitepress.render(
            io, mime, node, node.children, page, doc; kwargs...,
        )
    end
end

# ── Stopgap: heading anchors that contain LaTeX ──────────────────────────────
# DocumenterVitepress builds each heading as `## <text> {#<slug>}`, where the
# slug is Documenter's anchor label passed through its own
# `sanitized_anchor_label` — whose comment says "vitepress doesn't like special
# markdown characters in the id slug", but which only strips `[ ] ( ) *`.
#
# A heading whose text carries LaTeX yields a slug holding a backslash and
# braces. VitePress's `{#...}` parser rejects them, so it treats the whole suffix
# as *text*: the heading renders with the raw `{#...}` visible, the formula is
# dropped, and the same garbage lands in the "On this page" outline.
#
# Stripping those characters from the slug is safe: this narrows to headings
# only, leaving docstring anchors — which legitimately carry braces, are emitted
# as raw `<a id=…>`, and *are* linked to — untouched.
#
# Remove once `sanitized_anchor_label` covers these characters upstream.
function DocumenterVitepress.render(
        io::IO,
        mime::MIME"text/plain",
        node::Documenter.MarkdownAST.Node,
        header::Documenter.AnchoredHeader,
        page,
        doc;
        kwargs...,
    )
    anchor = header.anchor
    label = DocumenterVitepress.sanitized_anchor_label(anchor)
    id = replace(replace(label, r"[\\{}]" => ""), " " => "-")
    heading = first(node.children)
    println(io)
    print(io, "#"^(heading.element.level), " ")
    heading_iob = IOBuffer()
    DocumenterVitepress.render(heading_iob, mime, node, heading.children, page, doc; kwargs...)
    print(io, rstrip(String(take!(heading_iob))))
    print(io, " {#$(id)}")
    if haskey(kwargs, :inventory)
        item = DocumenterVitepress.InventoryItem(
            name = anchor.id,
            domain = "std",
            role = "label",
            dispname = DocumenterVitepress._get_inventory_dispname(
                anchor.id, Documenter.MDFlatten.mdflatten(anchor.node)
            ),
            priority = -1,
            uri = DocumenterVitepress._get_inventory_uri(doc, page, id),
        )
        push!(kwargs[:inventory], item)
    end
    println(io)
    return nothing
end

# ── Stopgap: ordered lists start at 2, and swallow their first item ──────────
# DocumenterVitepress numbers ordered-list items with `bullet(i) = "$(i+1). "`,
# but `enumerate` is already 1-based, so every ordered list comes out numbered
# from 2. It also emits no blank line before the list.
#
# Together those two do real damage: a list whose first marker is `2.` cannot
# interrupt a paragraph — CommonMark allows that only for a list starting at
# `1.` — so with no separating blank line the first item is absorbed into the
# preceding prose as plain text and the list begins at `3.`.
#
# Remove once the numbering is fixed upstream.
function DocumenterVitepress.render(
        io::IO,
        mime::MIME"text/plain",
        node::Documenter.MarkdownAST.Node,
        list::Documenter.MarkdownAST.List,
        page,
        doc;
        kwargs...,
    )
    bullet(i) = list.type === :ordered ? "$(i). " : "- "
    println(io)
    iob = IOBuffer()
    for (i, item) in enumerate(node.children)
        DocumenterVitepress.render(
            iob, mime, item, item.children, page, doc; prenewline = false, kwargs...
        )
        eachline = split(String(take!(iob)), '\n')
        # Continuation lines must line up under the marker's full width. Upstream
        # hard-codes two spaces, which fits `- ` but not `1. `: a display equation
        # inside an ordered item falls out of the list, splitting it in two.
        pad = " "^length(bullet(i))
        eachline[2:end] .= pad .* eachline[2:end]
        final_string = join(eachline, '\n')
        endswith(final_string, '\n') || (final_string *= "\n")
        print(io, bullet(i))
        print(io, final_string)
    end
    return nothing
end

makedocs(;
    # `clean = false` let pages deleted from the source survive in `build/` and
    # go on being deployed. Nothing writes into `build/` before `makedocs`, so
    # wiping it costs nothing.
    clean = true,
    modules = [DECUHR],
    remotes = nothing,
    authors = "Jean-François Barthélémy",
    sitename = "DECUHR.jl",
    # The favicon and the logo are picked up automatically from `docs/src/assets`,
    # and the sidebar is derived from `pages` below, so neither needs declaring.
    # The former `assets = ["assets/custom.css"]` is gone too: that stylesheet is
    # ported to `docs/src/.vitepress/theme/overrides.css`, which the theme loads
    # after its own `style.css`.
    format = DocumenterVitepress.MarkdownVitepress(;
        repo = "https://github.com/MicroPoroChemoMechanics/DECUHR.jl",
        devbranch = "main",
        devurl = "dev",
        deploy_url = "https://MicroPoroChemoMechanics.github.io/DECUHR.jl",
        description = "Adaptive cubature for vertex singularities in Julia",
    ),
    plugins = [bib],
    pages = [
        "Home" => "index.md",
        "Algorithm" => "algorithm.md",
        "Examples" => "generated/examples.md",
        "API Reference" => "api.md",
        "License" => "license.md",

        "References" => "references.md",
    ],
    checkdocs = :exports,
    # No blanket `warnonly`: with `checkdocs = :exports` the missing-docs check
    # is already narrow, and every exported name is on a page. Leaving the
    # exemption in place would only hide the next cross-reference that breaks —
    # and under VitePress a broken reference becomes a dead link that fails the
    # build with a Rollup stack trace instead of a named file and line.
)

# DocumenterVitepress writes a real directory per version rather than the
# symlinks Documenter used, so it needs its own `deploydocs`.
DocumenterVitepress.deploydocs(;
    repo = "github.com/MicroPoroChemoMechanics/DECUHR.jl.git",
    target = joinpath(@__DIR__, "build"),
    branch = "gh-pages",
    devbranch = "main",
    push_preview = false,
)
