using StaticArrays

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

sk!(out, H::SKHamiltonian{4}, U, bonds) = OldSK._sk4!(U, bonds, out)
# sk_d!{IO}(H::SKHamiltonian{IO, 4}, r, R, b, db, dout) = _sk4_d!(R/r, r, b, db, dout)

######## spd-orbital model

sk!(out, H::SKHamiltonian{9}, U, bonds) = OldSK._sk9!(U, bonds, out)
# sk_d!{IO}(H::SKHamiltonian{IO, 9}, r, R, b, db, dout) = _sk9_d!(R/r, r, b, db, dout)



########################################################################
