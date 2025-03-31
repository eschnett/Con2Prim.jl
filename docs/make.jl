# Generate documentation with this command:
# (cd docs && julia make.jl)

push!(LOAD_PATH, "..")

using Documenter
using Con2Prim

makedocs(; authors="Erik Schnetter", format=Documenter.HTML(), sitename="Con2Prim")

deploydocs(; devbranch="main", push_preview=true, repo="github.com/eschnett/Con2Prim.jl.git")
