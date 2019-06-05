
module NewAssembly


"""
Store a TB Hamiltonian in RI format (i.e. the bond integrals)
"""
mutable struct OffsiteHamiltonianRI{T}
   H::SHK
   Iat::Vector{Int}            # first atom index i
   Jat::Vector{Int}            # second atom index j
   R::Vector{SVector{3,T}}     # X[j] - X[i]
   V::Matrix{T}                # a bond integral value
   Nat::Int
end


Base.length(H::OffsiteHamiltonianRI) = length(H.V)

Base.iterate(H::OffsiteHamiltonianRI) = zip(H.Iat, H.Jat, H.R, H.V, H.bond_type)

# TODO:
# - k-dependence
# - overlaps


function Base.Matrix(H::OffsiteHamiltonianRI{ORB, T}) where {ORB, T}
   # total number of orbitals = size of matrix
   norb_local = length(H.orbitals)
   norb_global = H.Nat * norb_local
   Hout = zeros(T, norb_global, norb_global)
   Hloc = zeros(T, norb_local, norb_loca)
   # loop over all bond-integrals between any two atoms
   for (iat, jat, R, V, bond_type) in H
      @assert jat > iat # assume this convention for now
      # assemble the relevant sk block and the associated _local_ matrix Hloc
      # (this includes the symmetries within the local block!)
      # sig ∈ {1, -1} is the sign of the transposed element
      irow, jcol, sig = sk_block!(Hloc, R, bond_type, H.orbitals)
      # get the global indices from the local and the atom indices
      irowg = global_indices(H.orbitals, iat, irow)
      jcolg = global_indices(H.orbitals, jat, jcol)
      # write information into global H
      for (ig, jg, i, j) in zip(irowg, jcolg, irow, jcol)
         H[ig, jg] += V * Hloc[i, j]
         H[jg, ig] += sig * V * Hloc[i, j]
      end
   end
   return H
end


end
