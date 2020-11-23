
using StaticArrays
using JuLIP: Atoms, neighbourlist, cutoff
using NeighbourLists: max_neigs, sites
using LinearAlgebra: norm

"""
Store a TB Hamiltonian in RI format (i.e. the bond integrals)
"""
mutable struct SKHamiltonianRI{T}
   skh::SKH
   model::SKModel
   Nat::Int
   # ----------- off-diagonal blocks
   Iat::Vector{Int}            # first atom index i
   Jat::Vector{Int}            # second atom index j
   R::Vector{SVector{3,T}}     # X[j] - X[i]
   V::Matrix{T}                # bond integral values
   # ----------- diagonal blocks (require they  are diagonal for now!)
   D::Matrix{T}                # diagonal entries Nat x norbitals
end

Base.length(H::SKHamiltonianRI) = length(H.Iat)


# -----
# Prototypes

function eval_bonds! end
function eval_diag! end
function filter_nhd end
function getskh end

# ----------------------------------------------------------------
#    Two-Centre Case

function assembleRI(at::Atoms, model::TwoCentreModel)
   nlist = neighbourlist(at, cutoff(model))
   skh = getskh(model)

   # copy / allocate space
   Nat = length(at)
   Iat = copy(nlist.i)
   Jat = copy(nlist.j)
   # R = copy(nlist.R)
   R = [ ijR[3] for ijR in pairs(nlist) ]
   V = zeros(length(R), nbonds(skh))

   # call the function barriers
   temp = MVector(zeros(nbonds(skh))...)
   _assembleRI!(V, R, model, temp)

   # and for the diagonal entries ...
   D = zeros(Nat, nbonds(skh))
   _assembleRI_diag!(D, model, nlist, temp)

   return SKHamiltonianRI(skh, model, Nat, Iat, Jat, R, V, D)
end

function _assembleRI!(V, R, model::TwoCentreModel, temp)
   for n = 1:length(R)
      eval_bonds!(temp, model, norm(R[n]))
      V[n, :] .= temp
   end
   return V
end

# TODO: at the moment, this is generic -> next allow SKmodel to have
#       general diagonal blocks
function _assembleRI_diag!(D, model::SKModel, nlist, temp)
   for (i, j, r, R) in sites(nlist)
      eval_diag!(temp, model, r, R)
      D[i, :] .= temp
   end
   return D
end




# ----------------------------------------------------------------
#    General Case

function assembleRI(at::Atoms, model::SKModel)
   nlist = neighbourlist(at, offsite_cutoff(model))
   skh = getskh(model)

   # copy / allocate space
   Nat = length(at)
   nb = nbonds(skh)
   Iat = zeros(Int, nb)
   Jat = zeros(Int, nb)
   R = zeros(JVecF, nb)
   V = zeros(length(R), nb)

   # call the function barrier
   temp = MVector(zeros(nbonds(skh))...)
   _assembleRI!(Iat, Jat, R, V, nlist, model, temp)

   # and for the diagonal entries ...
   D = zeros(Nat, nbonds(skh))
   _assembleRI_diag!(D, model, nlist, temp)

   return SKHamiltonianRI(skh, model, Nat, Iat, Jat, R, V, D)
end


function _assembleRI!(Iat, Jat, R, V, nlist, model, temp)
   idx = 0
   R_filtered = sizehint!(JVecF[], max_neigs(nlist))  # allocate without growing the array
   for (i, J, r, R) in sites(nlist)
      for (nj, j) in enumerate(J)
         Rji  = R[j]  # R[j] = X[j] - X[i]
         # filter -> Ks, Rs = [X[k]-X[i] for k in ...]
         for k in Jh
            Rki = R[k]
            if (k != j) && filter_nhd(Rji, Rki)
               push!(R_filtered, Rki)
            end
         end
         # evaluate the bonds
         fill!(temp, 0.0)
         eval_bonds!(temp, model, Rji, R_filtered)
         empty!(R_filtered)    # this keeps the memory allocated
         # write into V
         idx += 1
         V[idx, :] .= temp
      end
   end
   return V
end



# ----------------------------------------------------------------
#    Assemble the Hamiltonian Matrix in Cartesian format

"""
number of local orbitals
"""
get_norbloc(H::SKH) = sum( 2*get_l(o)+1 for o in H.orbitals )

"""
construct global Hamiltonian matrix indices from atom indices and
local orbital indices
"""
global_index(iat, iloc, nat, norbloc) = (iat-1) * norbloc .+ iloc

function local_indices(nb::Integer, skh::SKH)
   io1, io2 = skh.b2o[nb] # orbital indices
   I1, I2 = skh.locorbidx[io1], skh.locorbidx[io2]   # local block indices
   # convert to SVectors
   return SVector(I1...), SVector(I2...)
end

function Base.Matrix(HRI::SKHamiltonianRI{T}) where {T}
   # ------------------------------------------------------------------------
   # this "outer" matrix assembly is intentionally type unstable; it prepares
   # everything needed for a fast assembly in _assemble_full_inner!, which is the
   # function barrier beyond which type stability should be strictly enforced
   # ------------------------------------------------------------------------
   norbloc = get_norbloc(HRI.skh)
   norbglob = HRI.Nat * norbloc

   # allocate the global hamiltonian matrix
   H = zeros(norbglob, norbglob)

   for nb = 1:length(HRI.skh.bonds)
      # compute the local orbital indices corresponding to the current bond type
      # TODO make sure these indices are SVectors
      irowloc, jcolloc = local_indices(nb, HRI.skh)
      @assert ( (irowloc isa SVector{N, <: Integer} where {N}) &&
                (jcolloc isa SVector{N, <: Integer} where {N}) )
      # call the function barrier
      _assemble_full_inner!(H, HRI, nb, HRI.skh.bonds[nb], irowloc, jcolloc, norbloc)
   end

   _assemble_full_inner_diag!(H, HRI, SVector((1:norbloc)...))  # TODO: split into TwoCentre / General

   return H
end

# TODO: some ideas to speed this up
#    - irowloc, jcolloc could become type parameters
#    - the global indices could be incremented rather than recomputed
#    - compute only upper triangular part, copy the lower triangular part
function _assemble_full_inner!(H, HRI, nb, b, irowloc, jcolloc, norbloc)
   sig = sksign(b)
   sigt = sksignt(b)
   for (n, (iat, jat, R)) in enumerate( zip(HRI.Iat, HRI.Jat, HRI.R) )
      V = HRI.V[n, nb]
      U = R / norm(R)
      φ, θ = carttospher(U[1], U[2], U[3])
      E12 = CodeGeneration.sk_gen(b, φ, θ)
      irowglob = global_index(iat, irowloc, HRI.Nat, norbloc)
      jcolglob = global_index(jat, jcolloc, HRI.Nat, norbloc)
      H[irowglob, jcolglob] += sig * V * E12
      if irowloc != jcolloc # non-diagonal blocks must be copied
         irowglobt = global_index(iat, jcolloc, HRI.Nat, norbloc)
         jcolglobt = global_index(jat, irowloc, HRI.Nat, norbloc)
         H[irowglobt, jcolglobt] += sigt * V * E12'
      end
   end
   return H
end


function _assemble_full_inner_diag!(H, HRI, iloc)
   norbloc = length(iloc)
   for nat = 1:HRI.Nat
      iglob = global_index(nat, iloc, HRI.Nat, norbloc)
      for j = 1:norbloc
         H[iglob[j], iglob[j]] = HRI.D[nat, j]
      end
   end
   return H
end

# TODO:
# - k-dependence
# - overlaps
