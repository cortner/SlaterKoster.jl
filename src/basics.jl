using StaticArrays

# part of SlaterKoster.jl

_int_to_bond = Dict( 0 => "σ", 1 => "π", 2 => "δ", 3 => "φ", 5 => "γ" )
_bond_to_int = Dict( "σ" => 0, "π" => 1, "δ" => 2, "φ" => 3, "γ" => 5 )
_int_to_orb = Dict( 0 => "s", 1 => "p", 2 => "d", 3 => "f", 4 => "g" )
_orb_to_int = Dict( "s" => 0, "p" => 1, "d" => 2, "f" => 3, "g" => 4 )

orb_to_int(s::Union{Char, AbstractString}) = _orb_to_int[String(s)]
bond_to_int(s::Union{Char, AbstractString}) = _bond_to_int[String(s)]
int_to_orb(i::Integer) = _int_to_orb[i]
int_to_bond(i::Integer) = _int_to_bond[i]
allowed_orbitals() = collect(keys(_orb_to_int))
allowed_bonds() = collect(keys(_bond_to_int))

"""
maximal symmetry symbol index given two angular momentum numbers.
"""
max_symbol(l1,l2) = min(l1, l2)


"""
`bondintegral_index` :

 1 ssσ  -> (0, 0, 0) = (l1, l2, M) -> 1
 2 psσ
 3 ppσ
 4 ppπ
 5 dsσ
 6 dpσ
 7 dpπ
 8 ddσ
 9 ddπ
10 ddδ

TODO: there must be a cleverer implementation of this!
"""
function bondintegral_index(l1, l2, sym)
   if l1 < l2
      l1, l2 = l2, l1
   end
   @assert sym <= max_symbol(l1, l2)

   idx = 0
   for _l1 = 0:l1, _l2 = 0:_l1, _sym = 0:max_symbol(_l1, _l2)
      idx += 1
      if (_l1, _l2, _sym) == (l1, l2, sym)
         return idx
      end
   end
   @error("we shouldn't be here!")
end

"""
`orbital_index(l::Integer, m::Integer)` : old indexing implementation for
compatibility! Should move to `OrbitalIndices`.

     l  m
s    0  0    -> 1
pz   1  -1   -> 2
py   1  0    -> 3
px   1   1   -> 4
...

TODO: there must be a cleverer implementation of this!
"""
function orbital_index(l::Integer, m::Integer)
   @assert abs(m) <= l
   idx = 0
   for _l = 0:l, _m = -_l:_l
      idx += 1
      if (_l, _m) == (l, m)
         return idx
      end
   end
   @error("we shouldn't be here!")
end

"""
`Orbitals(orbs)` : returns a struct that provides local indexing information
for the list of orbitals provided, as well as the associated bond integrals.

### Input

`orbs` must be an iterable containing `AbstractString`s, where each string
specifies an orbital. The order in the vector provides the ordering of the
orbitals that the TB model will use to assemble the Hamiltonians.

### Output

An `Orbitals` object which stores the information which orbitals are present
on some atom, including precomputed information to enable efficient indexing
and looping over the orbitals.
"""
struct Orbitals{ORB}
   v::Val{ORB}
   str::String
   orbidx::Matrix{Int}
   bondidx::Array{3, Int}
   orbidx_str::Dict{String, Int} # TODO: convert to NamedTuple for faster access
end

Base.show(io, ::Orbitals{ORB}) where ORB =
      show(io, "Orbitals{$(string(ORB))}")

# make this a get-index??

orbital_index(oi::Orbitals, s::AbstractString, l, m) =
   orbital_index(oi, oi.orbidx[s], l, m)

orbital_index(oi::Orbitals, n::Integer, l, m) =
   oi.orbidx[n, 1+l+m]


function Orbitals(orbs::AbstractVector{<: AbstractString})
   orbs = collect(orbs)
   norbs = length(orbs)
   # check correctness of the list of orbitals `orbs`
   allorbs = allowed_orbitals()
   for s in orbs
      @assert 1 <= length(s) <= 2
      @assert string(s[end]) in allorbs
      if length(s) == 2
         @assert Int('0') <= Int(s[1]) <= Int('9')
      end
   end
   # maximum l value
   L = maximum( orb_to_int(s[end]) for s in orbs )
   # allocate an array of indices
   orbidx_str = Dict{String, Int}()
   orbidx = zeros(Int, norbs, 2*L+1)
   bondidx = zeros(Int, norbs, norbs, max_symbol(L, L))
   orbidx_ctr = 0
   for (is, s) in enumerate(orbs)
      orbidx_str[s] = is
      l = orb_to_int(s[end])
      for m = -l:l
         orbidx_ctr += 1
         orbidx[is, 1+m+l] = orbidx_ctr
      end
   end
   valorbs = Val(Symbol(prod(orbs)))
   return Orbitals(valorbs, orbidx, orbidx_str)
end


"""
A `Bond` or Bond Integral basically represents a Vssσ, V_spσ, V_ppπ, etc
object. It is primarily used for indexing operations.
"""
struct Bond
   orb1::Symbol   # first orbital (s, p, d,...) or (1s, 2s, ...)
   orb2::Symbol   # second orbital
   sym::Symbol    # symmetry      (σ, π, ...)
end

"""
A variant of `Bond` storing integers instead of symbols
for faster access. Needs an `Orbitals` object to convert from
`IntBond` to `Bond`
"""
struct IntBond
   orb1::Int8
   orb2::Int8
   sym::Int8
end

"""
Given a bond `b::Bond` or `b::IntBond` we need to obtain `l1, l2, sym`
in order to compute the relevant matrix block, as well as the local
row and column indices `irow, jcol` corresponding to that block.
"""
function bond2orb(b::IntBond, orbitals::Orbitals)
   
end
