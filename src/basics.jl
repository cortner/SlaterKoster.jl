using StaticArrays

export @sko_str, @skb_str, SKOrbital, SKBond, get_l, bond_to_idx

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
get_l(s::Union{Symbol, Char}) = _orb_to_L[Symbol(s)]
get_l(s::AbstractString) = _orb_to_L[Symbol(s[end])]

for (sym, val) in _orb_to_L
   "get_l(::Val{:$sym}) = $val" |> Meta.parse |> eval
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
only the last letter will be used in determining the type of orbital
"""
struct SKOrbital{S}
   valS::Val{S}     # symbol of the orbital
   l::Int           # l-value
   str::String      # description
   idx::Int         # orbital index
end

Base.String(o::SKOrbital) = "sk:" * o.str

Base.show(io::IO, o::SKOrbital) = write(io, String(o))

function SKOrbital(str, idx::Integer = 0)
   @assert 1 <= length(str) <= 2
   @assert Symbol(str[end]) in allowed_orbitals()
   return SKOrbital(Val(Symbol(str[end])), get_l(str), str, idx)
end

SKOrbital(o::SKOrbital, idx::Integer) = SKOrbital(o.valS, o.l, o.str, idx) 

macro sko_str(str) SKOrbital(str) end

get_l(o::SKOrbital) = o.l # get_l(o.valS)
get_idx(o::SKOrbital) = o.idx
bondtypes(o1::SKOrbital, o2::SKOrbital) = bondtypes(get_l(o1), get_l(o2))

# import Base.==
# ==(o1::SKOrbital, o2::SKOrbital) = (o1.valS == o2.valS) && (o1.valN == o2.valN)

import Base: isless
isless(o1::SKOrbital, o2::SKOrbital) = (
      (get_l(o1), get_idx(o1)) < (get_l(o2), get_idx(o2)) )

# ----------------------------------------------------------------------------
#    Bond Implementation

struct SKBond{O1,O2,SYM}
   o1::SKOrbital{O1}
   o2::SKOrbital{O2}
   valSYM::Val{SYM}
end

Base.String(b::SKBond{O1,O2,SYM}) where {O1,O2,SYM} = "sk:$O1$O2$SYM"

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
   return SKBond(SKOrbital(String(o1)), SKOrbital(String(o2)), Val(sym))
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

#FHISigns
function sksignmat(l1, l2)
   l1_len = 2 * l1 + 1
   l2_len = 2 * l2 + 1
   M = ones(l1_len, l2_len)
   for i=1:l1_len
      for j=1:l2_len
         M[i,j] *= sksign(i - l1 - 1, 0)
         M[i,j] *= sksign(j - l2 - 1, 0)
      end
   end
   return M
end

sksignmatt(l1, l2) = sksignmat(l2, l1)
sksignmat(b::SKBond) = sksignmat(get_l(b)...)
sksignmatt(b::SKBond) = sksignmatt(get_l(b)...)

"""
`SKModel` : Most general Slater-Koster Tight-Binding model supertype.
"""
abstract type SKModel end

"""
`TwoCentreModel` : abstract subtype of `SKModel` from which all
2-centre models should be derived, for simplified assembly of
the hamiltonians.
"""
abstract type TwoCentreModel <: SKModel end
