module SlaterKoster

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

# define non-local types and functions
include("prototypes.jl")

# evaluate the SK blocks
include("skcore.jl")

# hamiltonian matrix assembly
include("matrixassembly.jl")


end # module
