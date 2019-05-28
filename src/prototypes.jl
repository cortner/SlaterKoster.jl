using StaticArrays

using JuLIP: AbstractAtoms, JVecF 

"""
`abstract type SKHamiltonian{NORB}` : supertype for Hamiltonian
operators of Slater-Koster type
"""
abstract type SKHamiltonian{NORB} end


# TODO: general implementation
"""
`nbonds(H)` : number of hopping functions needed to evaluate the
hamiltonian.
"""
nbonds(::SKHamiltonian{1}) = 1
nbonds(::SKHamiltonian{4}) = 4
nbonds(::SKHamiltonian{9}) = 10


"""
`hop(H::SKHamiltonian, r, i)`: where `r` is real, `i` integer, this
should return the
"""
hop

# TODO: we don't know yet what this means!
# hop_d(H::SKHamiltonian, r, i) = ForwardDiff.derivative(s -> hop(H,s,i), r)


"""
`hop!(H::SKHamiltonian, r, bonds)` : fill `bonds` with `hop(H, r, i)`
"""
function hop!(H::SKHamiltonian, r, bonds)
   for i = 1:nbonds(H)
      bonds[i] = hop(H, r, i)
   end
   return bonds
end

# # TODO: hop_d should return hop and hop_d, since we need both!!!!
# function hop_d!(H::SKHamiltonian, r, b, db)
#    for i = 1:nbonds(H)
#       b[i] = hop(H, r, i)
#       db[i] = hop_d(H, r, i)
#    end
#    return b, db
# end


"""
`overlap(H::SKHamiltonian, r, i)`: where `r` is real, `i` integer, this
should return the overlap SK parameter
"""
function overlap end

# overlap_d(H::SKHamiltonian, r::Real, i) = ForwardDiff.derivative(s->overlap(H,s,i), r)

"""
`overlap!(H::SKHamiltonian, r, bonds)` : fill `bonds` with `overlap(H, r, i)`
"""
function overlap!(H::SKHamiltonian, r, bonds)
   for i = 1:nbonds(H)
      bonds[i] = overlap(H, r, i)
   end
   return bonds
end

# function overlap_d!(H::SKHamiltonian, r, b, db)
#    for i = 1:nbonds(H)
#       b[i] = overlap(H, r, i)
#       db[i] = overlap_d(H, r, i)
#    end
#    return b, db
# end

# TODO: discuss the onsite terms
# rewrite onsite! as diagonals???
# but we don't provide AD for these since they depend on too many variables,
# so AD will necessarily be inefficient.
# function onsite! end
# function onsite_grad! end



# Matrix Assembly

const SKBlock{NORB} = SMatrix{NORB, NORB, Float64}

"""
`struct SparseSKH` :
a triplet sparse matrix kind of thing that stores pre-computed
hamiltonian blocks, and can efficiently generate k-dependent Hamiltonians,
either sparse of full.

### Methods:
* `full(A::SparseSKH, k)`: returns full `H`, `M` for given `k`-vector.
* `full!(out, A, k)`: same as `full`, but in-place

### Methods to be implemented:
* `sparse(A::SparseSKH, k)`: same but sparse
* `collect(A::SparseSKH, k)`: chooses full (small systems) or sparse (large
                           systems) based on some simple heuristic

### Fields:
* `H` : the hamiltonian used to construct it
* `at` : the associated atoms object
* `i, j` : row and column indices
* `first` : `first[n]` is the index in `i, j` for which `i[idx] == n`
* `vH` : Hamiltonian blocks
* `vM` : overlap blocks
* `Rcell` : each hamiltonian block is associated with an e^{i k ⋅ S} multiplier
            from the Bloch transform; this S is stored in Rcell.
"""
struct SparseSKH{HT, TV}  # v0.6: require that TV <: SKBlock{NORB}
   H::HT
   at::AbstractAtoms
   i::Vector{Int32}
   j::Vector{Int32}
   first::Vector{Int32}
   vH::Vector{TV}
   vM::Vector{TV}
   Rcell::Vector{JVecF}
end

Base.length(skh::SparseSKH) = length(skh.i)
