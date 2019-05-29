

using StaticArrays

using JuLIP: AbstractAtoms, positions, cutoff, neighbourlist
using NeighbourLists: sites, npairs
using LinearAlgebra: dot, diagind, Hermitian

export hamiltonian

"""
`indexblock`: auxiliary function to compute indices of Slater Koster orbitals,
this is returned as an SVector, i.e. it is generated on the stack so that no
heap memory is allocated.
"""
indexblock(n::Integer, H::SKHamiltonian{NORB}) where {NORB} =
   SVector{NORB, Int}( ((n-1)*NORB+1):(n*NORB) )

ndofs(H::SKHamiltonian{NORB}, at::AbstractAtoms) where {NORB} =
      ndofs(NORB, length(at))
ndofs(norb::Integer, nat::Integer) =
      norb * nat

function SparseSKH(H::SKHamiltonian{NORB}, at::AbstractAtoms) where {NORB}

   # here is a little hack that will turn an abstract type into a concrete type
   SKB = typeof(zero(SKBlock{NORB}))

   # get back positions since we still can't do without to get the cell-shifts
   # TODO: >>> maybe implement an alternative in JuLIP?
   X = positions(at)

   # ALLOCATIONS: TODO: move to a different method?
   # neighbourlist
   nlist = neighbourlist(at, cutoff(H))
   #      off-diagonal + diagonal
   nnz = npairs(nlist) + length(at)
   # allocate space for ordered triplet format
   i = zeros(Int32, nnz)
   j = zeros(Int32, nnz)
   first = zeros(Int32, length(at))
   vH = zeros(SKB, nnz)
   vM = zeros(SKB, nnz)
   Rcell = zeros(JVecF, nnz)
   # index into the triplet format
   idx = 0

   # allocate space to assemble the hamiltonian blocks, we use MMatrix
   # here, but is this really necessary? Matrix should do fine?
   H_nm = zero(MMatrix{NORB, NORB, Float64})
   M_nm = zero(MMatrix{NORB, NORB, Float64})
   bonds = zeros(nbonds(H))     # temporary array for storing the potentials

   # loop through sites
   for (n, neigs, r, R) in sites(nlist)
      first[n] = idx+1          # where do the triplet entries for atom n start?

      # --------- on-site terms ----------
      # add the diagonal/on-site entries
      onsite!(H, r, R, H_nm, M_nm)
      # onsite_overlap!(H, r, R, M_nm)   # TODO: should this be just (H, M)
      idx += 1
      i[idx], j[idx], vH[idx], vM[idx] = n, n, SKB(H_nm), SKB(M_nm)
      # Rvector corresponding to this block is just X[n]-X[n] = 0
      Rcell[idx] = zero(JVecF)

      # loop through the neighbours of the current atom (i.e. the bonds)
      for m = 1:length(neigs)
         U = R[m]/r[m]
         # hamiltonian block; TODO: new framework plugs into here!!!! <<<<<<<<<< !!!!!!
         sk!(H_nm, H, U, hop!(H, r[m], bonds))
         sk!(M_nm, H, U, overlap!(H, r[m], bonds))
         idx += 1
         i[idx], j[idx], vH[idx], vM[idx] = n, neigs[m], SKB(H_nm), SKB(M_nm)
         # compute the Rcell vector for these blocks
         Rcell[idx] = R[m] - (X[neigs[m]] - X[n])
      end
   end

   return SparseSKH(H, at, i, j, first, vH, vM, Rcell)
end


_alloc_full(skh::SparseSKH) = _alloc_full(skh.H, skh.at)

_alloc_full(H::SKHamiltonian, at::AbstractAtoms) =
         ( Matrix{ComplexF64}(undef, ndofs(H, at), ndofs(H, at)),
           Matrix{ComplexF64}(undef, ndofs(H, at), ndofs(H, at)) )

Base.Matrix(H::SparseSKH, k::AbstractVector = zero(JVecF)) =
         _full!(_alloc_full(H), H, k)

_full!(out::Tuple, H::SparseSKH, k::AbstractVector = zero(JVecF)) =
         _full!(out[1], out[2], H, k, H.H)

function _full!(Hout, Mout, skh, k, H::SKHamiltonian{NORB}) where {NORB}
   fill!(Hout, 0.0)
   fill!(Mout, 0.0)
   k = JVecF(k)
   for (i, j, H_ij, M_ij, S) in zip(skh.i, skh.j, skh.vH, skh.vM, skh.Rcell)
      eikR = exp( im * dot(k, S) )
      Ii, Ij = indexblock(i, H), indexblock(j, H)
      @inbounds for a = 1:NORB, b = 1:NORB
         Hout[Ii[a], Ij[b]] += H_ij[a, b] * eikR
         Mout[Ii[a], Ij[b]] += M_ij[a, b] * eikR
      end
   end
   # TODO: enforce numerical self-adjointness here???
   #       e.g. by wrapping into Hermition???
   Hout[diagind(Hout)] = real(Hout[diagind(Hout)])
   Mout[diagind(Mout)] = real(Mout[diagind(Mout)])
   return Hout, Mout
end


hamiltonian(H::SKHamiltonian, at::AbstractAtoms, k = zero(JVecF)) =
            Hermitian.(Matrix(SparseSKH(H, at), k))
