using StaticArrays

using SlaterKoster.CodeGeneration: sk_gen!

########################################################################
######## Interface for computing SK Blocks     #########################
########################################################################

"""
`sk(H, args...)` : allocating version of `sk!`
"""
sk(H::SKHamiltonian{NORB}, args...) where {NORB} =
      sk!(zero(MMatrix{NORB,NORB,Float64}), H, args...)

"""
`sk!(out, H, U, bonds)` : main method for assembling SK blocks

TODO: write documentation
"""


# the following is just an interface and needs to be replaced of course

######## s-orbital model
sk!(out, H::SKHamiltonian{1}, U, bonds) = setindex!(out, bonds[1], 1)

# TODO: wait until needed
# function sk_d!{IO}(H::SKHamiltonian{IO, 1}, r, R, b, db, dH_nm)
#    # dH_nm is 3 x 1 x 1 so we can just index it linearly    (NORB = 1)
#    for a = 1:3
#       dH_nm[a] = db[1] * R[a] / r
#    end
#    return dH_nm
# end

######## sp-orbital model

sk!(out, H::SKHamiltonian{4}, U, bonds) = OldSK.sk4!(U, bonds, out)
# sk_d!{IO}(H::SKHamiltonian{IO, 4}, r, R, b, db, dout) = _sk4_d!(R/r, r, b, db, dout)

######## spd-orbital model

sk!(out, H::SKHamiltonian{9}, U, bonds) = OldSK._sk9!(U, bonds, out)
# sk_d!{IO}(H::SKHamiltonian{IO, 9}, r, R, b, db, dout) = _sk9_d!(R/r, r, b, db, dout)



########################################################################


norb2L(::Val{1}) = Val(0)    # s
norb2L(::Val{4}) = Val(1)    # p
norb2L(::Val{9}) = Val(2)    # d
norb2L(::Val{16}) =Val(3)    # f

sk!(out, H::SKHamiltonian{NORB}, U, bonds) where {NORB} =
            _sk!(out, norb2L(Val{NORB}()), U, bonds)

function _sk!(out, valL::Val{L}, U, V) where {L}
   φ, θ = carttospher(U[1], U[2], U[3])
   return sk_gen!(out, valL, V, φ, θ)
end


# this assumes that the coordinates are normalised
# TODO: α=φ, β=θ
"""
INPUTS: (x,y,z): This represents a position vector which is given using Cartesian coordinates.
RETURNS: The corresponding polar coordinates - radial, azimuthal, polar.
"""
function carttospher(x,y,z)
   β = arccos(z)
   if x != 0
      α = arctan2(y,x)
   else
      if y > 0
         α = π/2
      elseif y < 0
         α = - π/2
      else
         α = 0.0
      end
   end
   return α, β
end
