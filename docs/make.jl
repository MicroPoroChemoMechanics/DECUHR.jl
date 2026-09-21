using Documenter
using DECUHR
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
    pages = [
        "Home" => "index.md",
        "Algorithm" => "algorithm.md",
        "Examples" => "generated/examples.md",
        "API Reference" => "api.md",
        "License" => "license.md",
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
