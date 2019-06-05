using StaticArrays

export @sko_str, @skb_str

# part of SlaterKoster.jl

# this is just indexing the bond types
_bond_to_idx = Dict( :σ => 0, :π => 1, :δ => 2, :φ => 3, :γ => 5 )
# these are the l-values associated with an s, p,... orbital
_orb_to_L = Dict( :s => 0, :p => 1, :d => 2, :f => 3, :g => 4 )

_L_to_orb = Dict( [ val => key  for (key, val) in _orb_to_L ]... )
_idx_to_bond = Dict( [ val => key for (key, val) in _bond_to_idx]... )


orb_to_L(s::Union{Symbol, Char, AbstractString}) = _orb_to_L[Symbol(s)]
bond_to_idx(s::Union{Symbol, Char, AbstractString}) = _bond_to_idx[String(s)]
L_to_orb(i::Integer) = _L_to_orb[i]
idx_to_bond(i::Integer) = _idx_to_bond[i]
allowed_orbitals() = collect(keys(_orb_to_L))
allowed_bonds() = collect(keys(_bond_to_idx))

"""
`max_symbol_idx(l1,l2)` :
maximal bond symbol index given two angular momentum numbers.
"""
max_symbol_idx(l1,l2) = min(l1, l2)

"""
`bondtypes` : get a list of bond types (σ, π, etc) acting between
two bonds specified either by `SKBond` of an integer (`l`-value)
"""
function bondtypes(l1::Integer, l2::Integer)
   sym = Symbol[]
   for symidx = 0:max_symbol_idx(l1, l2)
      push!(sym, idx_to_bond(symidx))
   end
   return sym
end


"""
`get_l`: return the l value of an orbital. There are two variants:
* `get_l(::Union{Symbol, Char, AbstractString})` : uses `Dict` lookup
* `get_l(::Val{*})` where `*` is a symbol : uses static dispatch (zero overhead)
"""
get_l(s::Union{Symbol, Char, AbstractString}) = _orb_to_L[Symbol(s)]

for (sym, val) in _orb_to_L
   "get_l(::Val{:$sym}) = $val" |> Meta.parse |> eval
   # eval( :( get_l(::Val{$sym}) = $val ) )
end

"""
`get_bidx` : associate a bond-index to a bond type σ, π, ...
"""
function get_bidx end

for (sym, val) in _bond_to_idx
   "get_bidx(::Val{:$sym}) = $val" |> Meta.parse |> eval
end

# ----------------------------------------------------------------------------
#     Orbital Implementation

"""
`struct Orbital`:

An orbital is specified by a string or symbol consisting of numbers 1..9 and
letters `s,p,d,f,g`. The `_` character is ignored. A string specifying an
orbital may contain exactly one letter and at most one number, which must precede
the letter. E.g., the following are admissible descriptions:
```
   "s", "p", "1s", "2s", ...
```
"""
struct SKOrbital{S, N}
   valS::Val{S}   # symbol of the orbital
   valN::Val{N}   # shell index
end

Base.String(o::SKOrbital{S, N}) where {S, N} = replace("sk:$N$S", "0" => "")

Base.show(io::IO, o::SKOrbital) = write(io, String(o))

function SKOrbital(str)
   @assert 1 <= length(str) <= 2
   @assert Symbol(str[end]) in allowed_orbitals()
   if length(str) == 1
      n = 0
   else
      n = parse(Int64, str[1])
   end
   return SKOrbital(Val(Symbol(str[end])), Val(n))
end

macro sko_str(str) SKOrbital(str) end

get_l(o::SKOrbital) = get_l(o.valS)
get_idx(o::SKOrbital{S,N}) where {S,N} = N
bondtypes(o1::SKOrbital, o2::SKOrbital) = bondtypes(get_l(o1), get_l(o2))

# import Base.==
# ==(o1::SKOrbital, o2::SKOrbital) = (o1.valS == o2.valS) && (o1.valN == o2.valN)

import Base: isless
isless(o1::SKOrbital, o2::SKOrbital) = (
      (get_l(o1), get_idx(o1)) < (get_l(o2), get_idx(o2)) )

# ----------------------------------------------------------------------------
#    Bond Implementation

struct SKBond{O1,O2,SYM,N1,N2}
   o1::SKOrbital{O1, N1}
   o2::SKOrbital{O2, N2}
   valSYM::Val{SYM}
end

Base.String(b::SKBond{O1,O2,SYM,N1,N2}) where {O1,O2,SYM,N1,N2} =
      replace("sk:$N1$O1$N2$O2$SYM", "0" => "")

Base.show(io::IO, b::SKBond) = write(io, String(b))

macro skb_str(str) SKBond(str) end

function SKBond(str)
   @assert length(str) == 3 || length(str) == 5
   if length(str) == 3
      n1 = n2 = 0
      o1, o2 = Symbol(str[1]), Symbol(str[2])
      sym = Symbol(str[3])
   else
      n1, n2 = parse(Int64, str[1]), parse(Int64, str[3])
      o1, o2 = Symbol(str[2]), Symbol(str[4])
      sym = Symbol(str[5])
   end
   @assert o1 in allowed_orbitals()
   @assert o2 in allowed_orbitals()
   @assert sym in allowed_bonds()
   return SKBond(SKOrbital(Val(o1), Val(n1)),
                 SKOrbital(Val(o2), Val(n2)),
                 Val(sym))
end

function SKBond(o1::SKOrbital, o2::SKOrbital, sym::Symbol)
   @assert sym in allowed_bonds()
   return SKBond(o1, o2, Val(sym))
end


get_l(b::SKBond) = get_l(b.o1), get_l(b.o2)

get_bidx(b::SKBond) = get_bidx(b.valSYM)

import Base: isless
isless(b1::SKBond, b2::SKBond) = (
      (b1.o2, b1.o1, get_bidx(b1)) < (b2.o2, b2.o1, get_bidx(b2)) )


# ----------------------------
# Sign conventions

struct StandardSigns end
signmod(::Type{StandardSigns}, args...) = 1

struct FHISigns end
signmod(::Type{FHISigns}, l1, l2, m1, m2) =
   _codegen.signmatrix(l1, l2)[l1+m1+1, l2+m2+1]

# Standard SK Sign convention???
# TODO: Confirm this is ok as a sign convention!!!
sksign(l1, l2) = (isodd(l1+l2) && (l1 > l2)) ? -1 : 1
sksignt(l1, l2) = sksign(l2, l1)
sksign(b::SKBond) = sksign(get_l(b)...)
sksignt(b::SKBond) = sksignt(get_l(b)...)

# ----------------------------
#  OLD STUFF TO BE DELETED



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
   @assert sym <= max_symbol_idx(l1, l2)

   idx = 0
   for _l1 = 0:l1, _l2 = 0:_l1, _sym = 0:max_symbol_idx(_l1, _l2)
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

#
# """
# `Orbitals(orbs)` : returns a struct that provides local indexing information
# for the list of orbitals provided, as well as the associated bond integrals.
#
# ### Input
#
# `orbs` must be an iterable containing `AbstractString`s, where each string
# specifies an orbital. The order in the vector provides the ordering of the
# orbitals that the TB model will use to assemble the Hamiltonians.
#
# ### Output
#
# An `Orbitals` object which stores the information which orbitals are present
# on some atom, including precomputed information to enable efficient indexing
# and looping over the orbitals.
# """
# struct Orbitals{ORB}
#    v::Val{ORB}
#    str::String
#    orbidx::Matrix{Int}
#    bondidx::Array{3, Int}
#    orbidx_str::Dict{String, Int} # TODO: convert to NamedTuple for faster access
# end
#
# Base.show(io, ::Orbitals{ORB}) where ORB =
#       show(io, "Orbitals{$(string(ORB))}")
#
# # make this a get-index??
#
# orbital_index(oi::Orbitals, s::AbstractString, l, m) =
#    orbital_index(oi, oi.orbidx[s], l, m)
#
# orbital_index(oi::Orbitals, n::Integer, l, m) =
#    oi.orbidx[n, 1+l+m]
#
#
# function Orbitals(orbs::AbstractVector{<: AbstractString})
#    orbs = collect(orbs)
#    norbs = length(orbs)
#    # check correctness of the list of orbitals `orbs`
#    allorbs = allowed_orbitals()
#    for s in orbs
#       @assert 1 <= length(s) <= 2
#       @assert string(s[end]) in allorbs
#       if length(s) == 2
#          @assert Int('0') <= Int(s[1]) <= Int('9')
#       end
#    end
#    # maximum l value
#    L = maximum( orb_to_L(s[end]) for s in orbs )
#    # allocate an array of indices
#    orbidx_str = Dict{String, Int}()
#    orbidx = zeros(Int, norbs, 2*L+1)
#    bondidx = zeros(Int, norbs, norbs, max_symbol(L, L))
#    orbidx_ctr = 0
#    for (is, s) in enumerate(orbs)
#       orbidx_str[s] = is
#       l = orb_to_L(s[end])
#       for m = -l:l
#          orbidx_ctr += 1
#          orbidx[is, 1+m+l] = orbidx_ctr
#       end
#    end
#    valorbs = Val(Symbol(prod(orbs)))
#    return Orbitals(valorbs, orbidx, orbidx_str)
# end
#
#
# # """
# # A `Bond` or Bond Integral basically represents a Vssσ, V_spσ, V_ppπ, etc
# # object. It is primarily used for indexing operations.
# # """
# # struct Bond
# #    orb1::Symbol   # first orbital (s, p, d,...) or (1s, 2s, ...)
# #    orb2::Symbol   # second orbital
# #    sym::Symbol    # symmetry      (σ, π, ...)
# # end
# #
# # """
# # A variant of `Bond` storing integers instead of symbols
# # for faster access. Needs an `Orbitals` object to convert from
# # `IntBond` to `Bond`
# # """
# # struct IntBond
# #    orb1::Int8
# #    orb2::Int8
# #    sym::Int8
# # end
#
# # """
# # Given a bond `b::Bond` or `b::IntBond` we need to obtain `l1, l2, sym`
# # in order to compute the relevant matrix block, as well as the local
# # row and column indices `irow, jcol` corresponding to that block.
# # """
# # function bond2orb(b::IntBond, orbitals::Orbitals)
# #
# # end


abstract type SKModel end
