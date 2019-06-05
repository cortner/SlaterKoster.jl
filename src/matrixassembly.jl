
using StaticArrays
using JuLIP: Atoms, neighbourlist
using NeighbourLists: max_neigs

"""
Store a TB Hamiltonian in RI format (i.e. the bond integrals)
"""
mutable struct OffsiteHamiltonianRI{T}
   H::SKH
   model::SKModel
   Iat::Vector{Int}            # first atom index i
   Jat::Vector{Int}            # second atom index j
   R::Vector{SVector{3,T}}     # X[j] - X[i]
   V::Matrix{T}                # a bond integral value
   Nat::Int
end

Base.length(H::OffsiteHamiltonianRI) = length(H.Iat)


# -----
# Prototypes

function eval_bonds! end
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
   R = copy(nlist.R)
   V = zeros(length(R), nbonds(skh))

   # call the function barrier
   temp = MVector(zeros(nbonds(skh))...)
   _assembleRI!(V, R, model, temp)

   return OffsiteHamiltonianRI(skh, model, Iat, Jat, R, V, Nat)
end

function _assembleRI!(V, R, model::TwoCentreModel, temp)
   for n = 1:length(R)
      eval_bonds!(temp, model, norm(R[n]))
      V[n, :] .= temp
   end
   return V
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

   return OffsiteHamiltonianRI(skh, model, Iat, Jat, R, V, Nat)
end


function _assembleRI!(Iat, Jat, R, V, nlist, model, temp)
   idx = 0
   R_filtered = sizehint!(JVecF[], max_neigs(nlist))  # allocate without growing the array
   for (i, J, r, R) in sites(nlist)
      for (nj, j) in enumerate(J)
         Rji  = R[j]-R[i]
         # filter -> Ks, Rs = [R[k]-R[i] for k in ...]
         for k in J
            Rki = R[k]-R[i]
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

get_norbloc(H::SKH) = sum( 2*get_l(o)+1 for o in H.orbitals )

global_index(iat, iloc, nat, norbloc) = (iat-1) * norbloc .+ iloc

function Base.Matrix(HRI::OffsiteHamiltonianRI{T}) where {T}
   norbloc = get_norbloc(HRI.skh)
   norbglob = HRI.Nat * norbloc

   # allocate the global hamiltonian matrix
   H = zeros(norbglob, norbglob)

   for nb = 1:length(H.skh.orbitals)
      # call the function barrier
      _assemble_full_inner!(H, HRI, nb, b, todo, norbloc) # TODO make sure these indices are SVectors
   end

   return H
end

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
      @. H[irowglob, jcolglob] = sig * V * E12
      # + the symmetrised variations....
   end
end



# TODO:
# - k-dependence
# - overlaps
