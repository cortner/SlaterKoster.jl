using StaticArrays
using SlaterKoster: get_l, get_bidx, bond_to_idx, sksignmat

export @skh_str, SKH


# this assumes that the coordinates are normalised
# TODO: α=φ, β=θ
"""
INPUTS: (x,y,z): This represents a position vector which is given using Cartesian coordinates.
RETURNS: The corresponding polar coordinates - radial, azimuthal, polar.
"""
function carttospher(x,y,z)
   r = sqrt(x*x+y*y+z*z)
   θ = acos(z/r)
   if x != 0
      φ = atan(y, x)
   else
      φ = sign(y) * π * 0.5
   end
   return φ, θ
end

struct SKH{SIGN}
   orbitals::Vector{SKOrbital}    # this should become a `species => orbitals` list
   bonds::Vector{SKBond}       # this should become a `(spec1, spec2) => bonds list`
   b2o::Vector{Tuple{Int,Int}}
   locorbidx::Vector{Vector{Int}}
   sig::Type{SIGN}
end

Base.String(H::SKH) = ("sk:" * prod((o.str*":") for o in H.orbitals))[1:end-1]

Base.show(io::IO, H::SKH) = write(io, String(H))

import Base.==
==(H1::SKH, H2::SKH) = all( getproperty(H1, p) == getproperty(H2, p)
                            for p in fieldnames(SKH) )

nbonds(H::SKH) = length(H.bonds)
norbitals(H::SKH) = length(H.orbitals)

"""
`allbonds(orbitals)` : compute all the `SKBond`s between any two
`SKOrbital`s in the list `orbitals` and return as `Vector`.
"""
function allbonds(orbitals::Vector{<: SKOrbital})
   norb = length(orbitals)
   bonds = SKBond[]
   idx = 0
   for i1 = 1:norb, i2 = 1:norb
      o1, o2 = orbitals[i1], orbitals[i2]
      for sym in bondtypes(o1, o2)
         idx += 1
         push!(bonds, SKBond(o1, o2, sym, idx))
      end
   end
   return bonds
end


function SKH(orbitals::AbstractVector{<: SKOrbital},
             bonds::AbstractVector{<: SKBond},
             sig = StandardSigns)
   # construct local orbital -> index mapping
   locorbidx = Vector{Int}[]
   idx = 0
   for orb in orbitals
      len = 2 * get_l(orb) + 1
      push!(locorbidx, Vector{Int}(idx .+ (1:len)))
      idx += len
   end
   # construct local bonds -> orbitals connectivity
   b2o = Tuple{Int, Int}[]
   for b in bonds
      i1 = findfirst(o -> o == b.o1, orbitals)
      i2 = findfirst(o -> o == b.o2, orbitals)
      # insist on an ordering convention
      push!(b2o, (i1, i2))
   end
   return SKH(orbitals, bonds, b2o, locorbidx, sig)
end

_create_idx(orbitals::AbstractVector{<: SKOrbital}) =
   [ SKOrbital(o, idx) for (idx, o) in enumerate(orbitals) ]

function SKH(orbitals::AbstractVector{<: SKOrbital}, sig=StandardSigns)
   orbs_indx = _create_idx(orbitals)
   return SKH(orbs_indx, allbonds(orbs_indx), sig)
end

# experimental constructor for a TB model
function SKH(str::AbstractString, sig = StandardSigns)
   if '0' <= str[1] <= '9'
      @assert iseven(length(str))
      orbitals = [ SKOrbital(str[n:n+1]) for n in 1:2:length(str) ]
   else
      orbitals = [ SKOrbital(str[n:n])   for n = 1:length(str) ]
   end
   return SKH(orbitals, sig)
end

macro skh_str(s) SKH(s) end


max_locidx(H::SKH) = maximum(maximum(I) for I in H.locorbidx)

alloc_block(H::SKH) = zeros(max_locidx(H::SKH), max_locidx(H::SKH))

"""
`sk2cart:` assemble a Slater-Koster matrix block.

**Warning:** this is type-unstable and should not be used to assemble large
Hamiltonians.
"""
function sk2cart_other(H::SKH, R, V)
   φ, θ = carttospher(R[1], R[2], R[3])
   E = alloc_block(H)
   for (b, Vb, (io1, io2)) in zip(H.bonds, V, H.b2o)
      G12 = CodeGeneration.sk_gen(b, φ, θ)
      I1 = H.locorbidx[io1]
      I2 = H.locorbidx[io2]
      E[I1, I2] .+= (sksign(b) * Vb) * G12
   end
   return E
end

function sk2cart_FHIaims(H::SKH, R, V)
   φ, θ = carttospher(R[1], R[2], R[3])
   E = alloc_block(H)
   for (b, Vb, (io1, io2)) in zip(H.bonds, V, H.b2o)
      G12 = CodeGeneration.sk_gen(b, φ, θ)
      I1 = H.locorbidx[io1]
      I2 = H.locorbidx[io2]
      E[I1, I2] .+= (sksign(b) * Vb) * G12 .* sksignmat(b)
   end
   return E
end

#sk2cart(H::SKH, R, V) = sk2cart_other(H::SKH, R, V)
sk2cart(H::SKH, R, V; FHIaims=false) = FHIaims ? sk2cart_FHIaims(H::SKH, R, V) : sk2cart_other(H::SKH, R, V)

function sk2cart_onsite(H::SKH, Rlist, Vlist)
   E = alloc_block(H)
   for (R, V) in zip(Rlist, Vlist)
       φ, θ = carttospher(R[1], R[2], R[3])
       for (b, Vb, (io1, io2)) in zip(H.bonds, V, H.b2o)
           G12 = CodeGeneration.sk_gen(b, φ, θ)
           I1 = H.locorbidx[io1]
           I2 = H.locorbidx[io2]
           E[I1, I2] .+= (sksign(b) * Vb) * G12
       end
   end
   return E
end

function sk2cart_onsite_FHIaims(H::SKH, Rlist, Vlist)
   E = alloc_block(H)
   for (R, V) in zip(Rlist, Vlist)
       φ, θ = carttospher(R[1], R[2], R[3])
       for (b, Vb, (io1, io2)) in zip(H.bonds, V, H.b2o)
           G12 = CodeGeneration.sk_gen(b, φ, θ)
           I1 = H.locorbidx[io1]
           I2 = H.locorbidx[io2]
           E[I1, I2] .+= (sksign(b) * Vb) * G12 .* sksignmat(b)
       end
   end
   return E
end

#sk2cart_onsite(H::SKH, R, V) = sk2cart_onsite_other(H::SKH, R, V)
sk2cart_onsite(H::SKH, R, V; FHIaims=false) = FHIaims ? sk2cart_onsite_FHIaims(H::SKH, R, V) : sk2cart_onsite_other(H::SKH, R, V)

function sk2cart_num_other(H::SKH, R, V)
   φ, θ = carttospher(R[1], R[2], R[3])
   E = alloc_block(H)
   for (b, Vb, (io1, io2)) in zip(H.bonds, V, H.b2o)
      G12 = CodeGeneration.sk_num(b, φ, θ)
      I1 = H.locorbidx[io1]
      I2 = H.locorbidx[io2]
      E[I1, I2] .+= (sksign(b) * Vb) * G12
   end
   return E
end

function sk2cart_num_FHIaims(H::SKH, R, V)
   φ, θ = carttospher(R[1], R[2], R[3])
   E = alloc_block(H)
   for (b, Vb, (io1, io2)) in zip(H.bonds, V, H.b2o)
      G12 = CodeGeneration.sk_num(b, φ, θ)
      I1 = H.locorbidx[io1]
      I2 = H.locorbidx[io2]
      E[I1, I2] .+= (sksign(b) * Vb) * G12 .* sksignmat(b)
   end
   return E
end

#sk2cart_num(H::SKH, R, V) = sk2cart_num_other(H::SKH, R, V)
sk2cart_num(H::SKH, R, V; FHIaims=false) = FHIaims ? sk2cart_num_FHIaims(H::SKH, R, V) : sk2cart_num_other(H::SKH, R, V)

"""
todo doc
"""
function cart2sk_other(H::SKH, R, E::AbstractArray)
   φ, θ = carttospher(R[1], R[2], R[3])
   V = zeros(length(H.bonds))
   for (I, (b, (io1, io2))) in enumerate(zip(H.bonds, H.b2o))
      G12 = CodeGeneration.sk_gen(b, φ, θ)
      I1 = H.locorbidx[io1]
      I2 = H.locorbidx[io2]
      b_l = get_bidx(b) # bond symbol to L
      val = sum(sksign(b) * E[I1, I2] .* G12) 
      if b_l > 0 # for bonds other than 's' or l>0
            V[I] += 0.5 * val 
      else # for 's' bond or l=0
            V[I] += val 
      end
   end
   return V
end

function cart2sk_FHIaims(H::SKH, R, E::AbstractArray)
   φ, θ = carttospher(R[1], R[2], R[3])
   V = zeros(length(H.bonds))
   for (I, (b, (io1, io2))) in enumerate(zip(H.bonds, H.b2o))
      G12 = CodeGeneration.sk_gen(b, φ, θ)
      I1 = H.locorbidx[io1]
      I2 = H.locorbidx[io2]
      b_l = get_bidx(b) # bond symbol to L
      val = sum(sksign(b) * E[I1, I2] .* G12 .* sksignmat(b)) 
      if prnt
          println("I,J:",I1,I2,"  VAL: ",val," G12:",G12," sign1:",sksign(b)," sign2:",sksignmat(b))
      end
      if b_l > 0 # for bonds other than 's' or l>0
            V[I] += 0.5 * val 
      else # for 's' bond or l=0
            V[I] += val 
      end
   end
   return V
end

#cart2sk(H::SKH, R, V) = cart2sk_other(H::SKH, R, V)
cart2sk(H::SKH, R, V; FHIaims=false) = FHIaims ? cart2sk_FHIaims(H::SKH, R, V) : cart2sk_other(H::SKH, R, V)

function cart2sk_num_other(H::SKH, R, E::AbstractArray)
   φ, θ = carttospher(R[1], R[2], R[3])
   V = zeros(length(H.bonds))
   for (I, (b, (io1, io2))) in enumerate(zip(H.bonds, H.b2o))
      G12 = CodeGeneration.sk_num(b, φ, θ)
      I1 = H.locorbidx[io1]
      I2 = H.locorbidx[io2]
      b_l = get_bidx(b) # bond symbol to L
      val = sum(sksign(b) * E[I1, I2] .* G12) 
      if b_l > 0 # for bonds other than 's' or l>0
            V[I] += 0.5 * val 
      else # for 's' bond or l=0
            V[I] += val 
      end
   end
   return V
end

function cart2sk_num_FHIaims(H::SKH, R, E::AbstractArray)
   φ, θ = carttospher(R[1], R[2], R[3])
   V = zeros(length(H.bonds))
   for (I, (b, (io1, io2))) in enumerate(zip(H.bonds, H.b2o))
      G12 = CodeGeneration.sk_num(b, φ, θ)
      I1 = H.locorbidx[io1]
      I2 = H.locorbidx[io2]
      b_l = get_bidx(b) # bond symbol to L
      val = sum(sksign(b) * E[I1, I2] .* G12 .* sksignmat(b)) 
      if b_l > 0 # for bonds other than 's' or l>0
            V[I] += 0.5 * val 
      else # for 's' bond or l=0
            V[I] += val 
      end
   end
   return V
end

#cart2sk_num(H::SKH, R, V) = cart2sk_num_other(H::SKH, R, V)
cart2sk_num(H::SKH, R, V; FHIaims=false) = FHIaims ? cart2sk_num_FHIaims(H::SKH, R, V) : cart2sk_num_other(H::SKH, R, V)
