module SlaterKoster

using Reexport

# Functionality we need - by order of
#  (1) correct s, p, d orbitals
#      INPUT: U, hop -> OUTPUT: matrix
#  (2) AD derivatives
#        - w.r.t. U
#        - w.r.t. hop
#        all further derivatives are handles outside
#  (3) correct f and higher orbitals
#  (4) fast evaluation
#  (5) fast derivatives

# my old codes from the previous TB package
include("oldsk.jl")

include("basics.jl")

# define non-local types and functions
include("prototypes.jl")

# evaluate the SK blocks
include("code_generation.jl")
include("skcore.jl")

# hamiltonian matrix assembly
include("matrixassembly.jl")
# exports: hamiltonian

# Kwon model - primarily for testing...
include("kwon.jl")
@reexport using SlaterKoster.Kwon

end # module
