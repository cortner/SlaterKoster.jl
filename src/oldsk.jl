
module OldSK

using StaticArrays, JuLIP
using JuLIP: AbstractAtoms

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
`isorthogonal(H) -> Bool` : specify whether the tb model is orthogonal
or non-orthogonal (has non-trivial overlap). Default is `false`.
"""
isorthogonal(H::SKHamiltonian) = false


"""
`hop(H::SKHamiltonian, r, i)`: where `r` is real, `i` integer, this
should return the
"""
function hop end

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
function onsite! end
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


"""
Assemble a hamiltonian block for a Slater-Koster type
hamiltonian, with 4 orbitals (s, p).

### Parameters:
* `U` : R / |R|, orientation of the bond, must be a 3-vector
* `hop` : vector of hopping function values, must be a 4-vector
* `mat` : output matrix, must be at least 4 x 4
"""
function sk4!(U, hop, mat)
   l, m, n = U[1], U[2], U[3]
   # 4 orbitals are s, px, py, pz, these are the mat-indices
   # 4 bond types are : ssσ,spσ,ppσ,ppπ, these are the hop-indices
   #                     1   2   3   4
   mat[1,1] = hop[1]                            # E_ss = V_ssσ
   mat[2,2] = l*l * hop[3] + (1-l*l) * hop[4]   # E_xx = l² V_ppσ + (1-l²) V_ppπ
   mat[3,3] = m*m * hop[3] + (1-m*m) * hop[4]   # E_yy = m² V_ppσ + (1-m²) V_ppπ
   mat[4,4] = n*n * hop[3] + (1-n*n) * hop[4]   # E_zz = n² V_ppσ + (1-n²) V_ppπ
   mat[1,2] = l * hop[2]                        # E_sx = l V_spσ
   mat[1,3] = m * hop[2]                        # E_sy = m V_spσ
   mat[1,4] = n * hop[2]                        # E_sz = n V_spσ
   mat[2,1] = - mat[1,2]                        # E_xs = - E_sx
   mat[3,1] = - mat[1,3]                        # E_ys = - E_sy
   mat[4,1] = - mat[1,4]                        # E_zs = - E_sz
   mat[2,3] = l * m * (hop[3] - hop[4])         # E_xy = l m (V_ppσ - V_ppπ)
   mat[2,4] = l * n * (hop[3] - hop[4])         # E_xz = l n (V_ppσ - V_ppπ)
   mat[3,4] = m * n * (hop[3] - hop[4])         # E_yz = n m (V_ppσ - V_ppπ)
   mat[3,2] =  mat[2,3]                         # E_yx = E_xy
   mat[4,2] =  mat[2,4]                         # E_zx = E_xz
   mat[4,3] =  mat[3,4]                         # E_zy = E_yz
   return mat
end



"""
Assemble a hamiltonian block for a slater-koster type
hamiltonian, with 9 orbitals (s, p, d).

* `U` : orientation of the bond
* `hop` : hopping functions for bond type, 10-vector
* `mat` : output matrix 9 x 9
"""
function sk9!(U, hop, mat)
   # fill the [1:4, 1:4] block
   sk4!(U, hop, mat)
   # and then all the rest
   l, m, n = U[1], U[2], U[3]
   # ll, mm, nn = l*l, m*m, n*n    # TODO: for faster evaluation???
   # lm, ln, mn = l*m, l*n, m*n
   # rt3 = √3

   # sd
   mat[1,5] = √3 * l * m * hop[5]
   mat[1,6] = √3 * m * n * hop[5]
   mat[1,7] = √3 * l * n * hop[5]
   mat[1,8] = √3/2 * (l^2 - m^2) * hop[5]
   mat[1,9] = ( n^2 - (l^2 + m^2)/2 ) * hop[5]
   mat[5,1] = mat[1,5]
   mat[6,1] = mat[1,6]
   mat[7,1] = mat[1,7]
   mat[8,1] = mat[1,8]
   mat[9,1] = mat[1,9]

   # pd
   mat[2,5] = √3 * l * l * m * hop[6] + m * (1.0 - 2.0 * l^2) * hop[7]
   mat[2,6] = √3 * l * m * n * hop[6] - 2.0 * l * m * n * hop[7]
   mat[2,7] = √3 * l * l * n * hop[6] + n * (1.0 - 2.0 * l^2) * hop[7]
   mat[2,8] = √3/2 * l * (l^2 - m^2) * hop[6] + l * (1.0 - l^2 + m^2) * hop[7]
   mat[2,9] = l * (n^2 - (l^2 + m^2)/2) * hop[6] - √3 * l * n^2 * hop[7]
   mat[5,2] = - mat[2,5]
   mat[6,2] = - mat[2,6]
   mat[7,2] = - mat[2,7]
   mat[8,2] = - mat[2,8]
   mat[9,2] = - mat[2,9]
   mat[3,5] = √3 * l * m * m * hop[6] + l * (1.0 - 2.0 * m^2) * hop[7]
   mat[3,6] = √3 * m * m * n * hop[6] + n * (1.0 - 2.0 * m^2) * hop[7]
   mat[3,7] = √3 * l * m * n * hop[6] - 2.0 * l * m * n * hop[7]
   mat[3,8] = √3/2 * m * (l^2 - m^2) * hop[6] - m * (1.0 + l^2 - m^2) * hop[7]
   mat[3,9] = m * (n^2 - (l^2 + m^2)/2) * hop[6] - √3 * m * n^2 * hop[7]
   mat[5,3] = - mat[3,5]
   mat[6,3] = - mat[3,6]
   mat[7,3] = - mat[3,7]
   mat[8,3] = - mat[3,8]
   mat[9,3] = - mat[3,9]
   mat[4,5] = √3 * l * m * n * hop[6] - 2.0 * l * m * n * hop[7]
   mat[4,6] = √3 * m * n * n * hop[6] + m * (1.0 - 2.0 * n^2) * hop[7]
   mat[4,7] = √3 * l * n * n * hop[6] + l * (1.0 - 2.0 * n^2) * hop[7]
   mat[4,8] = √3/2 * n * (l^2 - m^2) * hop[6] - n * (l^2 - m^2) * hop[7]
   mat[4,9] = n * (n^2 - (l^2 + m^2)/2) * hop[6] + √3 * n * (l^2 + m^2) * hop[7]
   mat[5,4] = - mat[4,5]
   mat[6,4] = - mat[4,6]
   mat[7,4] = - mat[4,7]
   mat[8,4] = - mat[4,8]
   mat[9,4] = - mat[4,9]

   # dd
   mat[5,5] = 3.0 * l^2 * m^2 * hop[8] + (l^2 + m^2 - 4.0 * l^2 * m^2) * hop[9] +
               (n^2 + l^2 * m^2) * hop[10]
   mat[6,6] = 3.0 * m^2 * n^2 * hop[8] + (m^2 + n^2 - 4.0 * m^2 * n^2) * hop[9] +
               (l^2 + m^2 * n^2) * hop[10]
   mat[7,7] = 3.0 * l^2 * n^2 * hop[8] + (l^2 + n^2 - 4.0 * l^2 * n^2) * hop[9] +
               (m^2 + l^2 * n^2) * hop[10]
   mat[8,8] = 3.0/4 * (l^2 - m^2)^2 * hop[8] + (l^2 + m^2 - (l^2 - m^2)^2) * hop[9] +
               (n^2 + (l^2 - m^2)^2 /4 ) * hop[10]
   mat[9,9] = (n^2 - (l^2 + m^2) /2)^2 * hop[8] + 3.0 * n^2 * (l^2 + m^2) * hop[9] +
               3.0/4 * (l^2 + m^2)^2 * hop[10]
   mat[5,6] = 3.0 * l * m^2 * n * hop[8] + l * n * (1.0 - 4.0 * m^2) * hop[9] +
               l * n * (m^2 - 1.0) * hop[10]
   mat[5,7] = 3.0 * l^2 * m * n * hop[8] + m * n * (1.0 - 4.0 * l^2) * hop[9] +
               m * n * (l^2 - 1.0) * hop[10]
   mat[5,8] = 3.0/2 * l * m * (l^2 - m^2) * hop[8] + 2.0 * l * m * (m^2 - l^2) * hop[9] +
               1.0/2 * l * m * (l^2 - m^2) * hop[10]
   mat[5,9] = √3 * l * m * (n^2 - (l^2 + m^2)/2) * hop[8] - 2.0*√3 * l * m * n^2 * hop[9] +
               √3/2 * l * m * (1.0 + n^2) * hop[10]
   mat[6,7] = 3.0 * l * m * n^2 * hop[8] + l * m * (1.0 - 4.0 * n^2) * hop[9] +
               l * m * (n^2 - 1.0) * hop[10]
   mat[6,8] = 3.0/2 * m * n * (l^2 - m^2) * hop[8] -
               m * n * (1.0 + 2.0 * (l^2 - m^2)) * hop[9] +
               m * n * (1.0 + (l^2 - m^2) /2) * hop[10]
   mat[6,9] = √3 * m * n * (n^2 - (l^2 + m^2)/2) * hop[8] +
               √3 * m * n * (l^2 + m^2 - n^2) * hop[9] -
               √3/2 * m * n * (l^2 + m^2) * hop[10]
   mat[7,8] = 3.0/2 * l * n * (l^2 - m^2) * hop[8] +
               l * n * (1.0 - 2.0 * (l^2 - m^2)) * hop[9] -
               l * n * (1.0 - (l^2 - m^2) /2) * hop[10]
   mat[7,9] = √3 * l * n * (n^2 - (l^2 + m^2) /2) * hop[8] +
               √3 * l * n * (l^2 + m^2 - n^2) * hop[9] -
               √3/2 * l * n * (l^2 + m^2) * hop[10]
   mat[8,9] = √3/2 * (l^2 - m^2) * (n^2 - (l^2 + m^2) /2) * hop[8] +
               √3 * n^2 * (m^2 - l^2) * hop[9] +
               √3/4 * (1.0 + n^2) * (l^2 - m^2) * hop[10]
   mat[6,5] = mat[5,6]
   mat[7,5] = mat[5,7]
   mat[8,5] = mat[5,8]
   mat[9,5] = mat[5,9]
   mat[7,6] = mat[6,7]
   mat[8,6] = mat[6,8]
   mat[9,6] = mat[6,9]
   mat[8,7] = mat[7,8]
   mat[9,7] = mat[7,9]
   mat[9,8] = mat[8,9]
   return mat
end


function sk4_d!(U, r, hh, dhh, dmat)
   # U = ∇r, so it appears here both as the orientation and as ∇r
   l, m, n = U[1], U[2], U[3]
   ll, lm, ln, mm, mn, nn = l*l, l*m, l*n, m*m, m*n, n*n
   dl = ( (1.0-ll)/r ,     - lm/r ,     - ln/r )
   dm = (     - lm/r , (1.0-mm)/r ,     - mn/r )
   dn = (     - ln/r ,     - mn/r , (1.0-nn)/r )

   for d = 1:3
      dmat[d,1,1] = dhh[1] * U[d]
      dmat[d,2,2] = (ll * dhh[3] + (1 .- ll) * dhh[4]) * U[d] + 2*l * hh[3] * dl[d] - 2*l * hh[4] * dl[d]
      dmat[d,3,3] = (mm * dhh[3] + (1 .- mm) * dhh[4]) * U[d] + 2*m * hh[3] * dm[d] - 2*m * hh[4] * dm[d]
      dmat[d,4,4] = (nn * dhh[3] + (1 .- nn) * dhh[4]) * U[d] + 2*n * hh[3] * dn[d] - 2*n * hh[4] * dn[d]
      dmat[d,1,2] = l * dhh[2] * U[d] + hh[2] * dl[d]
      dmat[d,1,3] = m * dhh[2] * U[d] + hh[2] * dm[d]
      dmat[d,1,4] = n * dhh[2] * U[d] + hh[2] * dn[d]
      dmat[d,2,1] = - dmat[d,1,2]
      dmat[d,3,1] = - dmat[d,1,3]
      dmat[d,4,1] = - dmat[d,1,4]
      dmat[d,2,3] = lm * (dhh[3] - dhh[4]) * U[d] + (dl[d] * m + l * dm[d]) * (hh[3] - hh[4])
      dmat[d,2,4] = ln * (dhh[3] - dhh[4]) * U[d] + (dl[d] * n + l * dn[d]) * (hh[3] - hh[4])
      dmat[d,3,4] = mn * (dhh[3] - dhh[4]) * U[d] + (dm[d] * n + m * dn[d]) * (hh[3] - hh[4])
      dmat[d,3,2] =  dmat[d,2,3]
      dmat[d,4,2] =  dmat[d,2,4]
      dmat[d,4,3] =  dmat[d,3,4]
   end
   return dmat
end



sk!(out, H::SKHamiltonian{4}, U, bonds) = sk4!(U, bonds, out)
sk!(out, H::SKHamiltonian{9}, U, bonds) = sk4!(U, bonds, out)


using StaticArrays

using JuLIP: AbstractAtoms, positions, cutoff, neighbourlist
using NeighbourLists: sites, npairs
using LinearAlgebra: dot, diagind, Hermitian


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
   # for (n, neigs, r, R) in sites(nlist)
   for (n, neigs, R) in sites(nlist)
      r = norm.(R)
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






using Parameters

using JuLIP.Potentials: fcut, fcut_d, SplineCutoff,
   PairPotential, ZeroPairPotential, EAM, @analytic

import JuLIP: cutoff

# export KwonHamiltonian

"""
`KwonHamiltonian <: SKHamiltonian{4}`

Hamiltonian for an orthogonal sp TB model of Si developed by Kwon et al [1].

This implementation deviates  from [1] in how the cut-off is applied:
instead of "patching" a cubic spline between r1 and rcut, we simply multiply
with a quintic spline on the interval [0.5 (rcut + r0), rcut].

[1] I. Kwon, R. Biswas, C. Z. Wang, K. M. Ho and C. M. Soukoulis.
Transferable tight-binding models for silicon.
Phys Rev B 49 (11), 1994.
"""
@with_kw struct KwonHamiltonian <: SKHamiltonian{4}
   r0::Float64 = 2.360352   # Å
   Es::Float64 = -5.25       # eV
   Ep::Float64 = 1.2        # eV
   E0::Float64 = 8.7393204  # eV   8.7393204  is the original value; but it is just a constant
   # ------------------------------- Electronic Parameters
   # α   1     2    3     4
   #    ssσ   spσ  ppσ   ppπ
   hr0::NTuple{4, Float64} = (-2.038, 1.745, 2.75, -1.075)
   nc::NTuple{4, Float64} = (9.5, 8.5, 7.5, 7.5)
   rc::NTuple{4, Float64} = (3.4, 3.55, 3.7, 3.7)
   # ------------------------------- 2-body parameters
   m::Float64 = 6.8755
   mc::Float64 = 13.017
   dc::Float64 = 3.66995     # Å
   C::NTuple{4, Float64} = (2.1604385, -0.1384393, 5.8398423e-3, -8.0263577e-5)
   # --------------------------------
   r1::Float64 = 3.260176     # Å  (start of cut-off) >>> different r1 from Kwon paper
   rcut::Float64 = 4.16       # Å  (end of cut-off)
end

cutoff(H::KwonHamiltonian) = H.rcut

isorthogonal(H::KwonHamiltonian) = true

kwon_hop(H::KwonHamiltonian, r::Real, α) = ( H.hr0[α] * (H.r0 / r)^2 *
               exp( - 2 * (r/H.rc[α])^H.nc[α] + 2 * (H.r0/H.rc[α])^H.nc[α] ) )

hop(H::KwonHamiltonian, r::Real, α) = kwon_hop(H, r, α) * fcut(r, H.r1, H.rcut)

overlap(H::KwonHamiltonian, r, i) = 0.0

function onsite!(H::KwonHamiltonian, _r, _R, H_nn, M_nn)
   fill!(H_nn, 0.0)
   fill!(M_nn, 0.0)
   H_nn[1,1] = H.Es
   H_nn[2,2] = H_nn[3,3] = H_nn[4,4] = H.Ep
   M_nn[1,1] = M_nn[2,2] = M_nn[3,3] = M_nn[4,4] = 1.0
   return H_nn, M_nn
end


end
