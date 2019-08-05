module SlaterKoster

using Reexport

# Functionality we need - by order of
#  (2) AD derivatives
#        - w.r.t. U
#        - w.r.t. hop
#        all further derivatives are handles outside
#  (3) correct f and higher orbitals : implemented -> need to check!
#  (4) fast evaluation
#  (5) fast derivatives

# my old codes from the previous TB package
#include("oldsk.jl")

# basic SK machinery
include("basics.jl")

# evaluate the SK blocks
include("code_generation.jl")
include("skcore.jl")

# hamiltonian matrix assembly
include("matrixassembly.jl")
# exports: hamiltonian

# Some basic models - mostly for testing?
include("models.jl")
@reexport using SlaterKoster.Models
# This exports
#  - KwonTB

end # module
