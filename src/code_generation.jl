
__precompile__(false)
module CodeGeneration

import JSON, Calculus
using StaticArrays

using PyCall: pyimport
using JuLIP: save_json
using SlaterKoster: max_symbol_idx
using SlaterKoster: SKBond, get_l, get_bidx, StandardSigns, sksign, signmod

# load the python part of the code generation
# const _pysys = pyimport("sys")
# push!(_pysys."path", @__DIR__())
# const _codegen = pyimport("codegen")

sympy = pyimport("sympy")
sp_spin = pyimport("sympy.physics.quantum.spin")

const spsin = sympy.sin
const spcos = sympy.cos
const Rotation = sp_spin.Rotation

_lookupkey(l1, l2, m1, m2, sym) = "$l1,$l2,$m1,$m2,$sym"
_fname_sktable() = joinpath(@__DIR__(), "sktable.json")

"""
`sk_table(L::Integer)` : read or create a table of SK matrix element expressions
"""
function sk_table(L::Integer)
   filepath = _fname_sktable()
   try
      tbl = JSON.parsefile(filepath)
      if tbl["L"] >= L
         # return the existing table
         return tbl
      end
      @info("""The existing lookup table does not contain a high enough degree.
               I'm now going to create a new one which can take a little while.""")
   catch
      @info("""Reading SK lookup table was unsuccesful, I'm now going to create
               a new one. This could take a while (O(minutes))""")
   end

   # create a new lookup table
   norb = orbital_index(L, L)
   tbl = Dict{String, Any}("L" => L)
   for l1 = 0:L, l2=l1:L, m1=-l1:l1, m2=-l2:l2
      for sym = 0:max_symbol_idx(l1,l2)
         tbl[_lookupkey(l1,l2,m1,m2,sym)] = Gsym(l1,l2,m1,m2,sym)
      end
   end
   save_json(filepath, tbl)

   # return the new table
   return tbl
end

"""
INPUTS: l1: Angular quantum number of atom 1.
        l2: Angular quantum number of atom 2.
        m1: Magnetic quantum number of atom 1.
        m2: Magnetic quantum number of atom 2.
        sym: A float which represents the type of symmetry bond we are considering, e.g. M = 0
           for sigma, M = 1 for pi, M = 2 for delta, etc.
        alpha: Azithmuthal coordinate.
        beta: Polar coordinate.
RETURNS: This function returns the relevant coefficient that should be are used in writing the
         Slater-Koster transformations.
"""
function Gsym(l1, l2, m1, m2, sym)

   φ, θ = sympy.symbols("phi, theta")

   τ(m) = (m >= 0) ? 1 : 0

   A(m, φ) = ( (m == 0) ? sympy.sqrt(2)/2 :
               (-1)^m * (τ(m) * spcos(abs(m)*φ) + τ(-m) * spsin(abs(m)*φ)) )

   B(m, φ) = ( (m == 0) ? 0 :
               (-1)^m * (τ(-m) * spcos(abs(m)*φ) - τ(m) * spsin(abs(m)*φ) ) )

   rot(l, m, sym, θ) = Rotation.d(l, m, sym, θ).doit()

   S(l, m, sym, φ, θ) = A(m, φ) * ( (-1)^sym * rot(l, abs(m),  sym, θ)
                                             + rot(l, abs(m), -sym, θ) )

   T(l, m, sym, φ, θ) = B(m, φ) * ( (-1)^sym * rot(l, abs(m),  sym, θ)
                                             - rot(l, abs(m), -sym, θ) )

   if sym == 0
      expr = (2 * A(m1, φ) * rot(l1, abs(m1), 0, θ)
                * A(m2, φ) * rot(l2, abs(m2), 0, θ) )
   else
      expr = (   S(l1, m1, sym, φ, θ) * S(l2, m2, sym, φ, θ)
               + T(l1, m1, sym, φ, θ) * T(l2, m2, sym, φ, θ) )
   end

   expr = expr.simplify()

   # x, y, z = sympy.symbols("x y z")
   # expr = expr.simplify()
   # expr = expr.subs(spsin(φ) * spsin(θ), y)
   # expr = expr.subs(spcos(φ) * spsin(θ), x)
   # expr = expr.subs(spcos(θ), z)
   # # expr = expr.subs(spsin(φ)^2, 1.0 - spcos(φ)^2)
   # expr = expr.simplify()

   # convert the sympy expression into a string
   return py2jlcode(expr.__str__())
end


function py2jlcode(str)
   # replace some symbols to turn python into julia code
   str = replace(str, "**" => "^")
   # φ : Azithmuthal coordinate
   str = replace(str, "phi" => "φ")   # the python code uses alpha here
   # θ : Polar angle
   str = replace(str, "theta" => "θ")   # the python code uses alpha here
end



"""
Generate a standard Slater-Koster off-diagonal block with inputs being
the bond integrals `V` and spherical coordinates `φ, θ`.

TODO: document the sign convention.
"""
@generated function sk_gen!(g, ::Val{L}, V, φ, θ,
                            sgnconv::Type{SGN} = StandardSigns
                            ) where {L, SGN}
   # get the SK expressions table
   tbl = sk_table(L)
   code = Expr[]
   for l1 = 0:L, l2 = l1:L, m1 = -l1:l1, m2 = -l2:l2
      # matrix indices, skip the lower-triangular part
      I1 = orbital_index(l1, m1)
      I2 = orbital_index(l2, m2)
      if I2 < I1; continue; end

      # start assembling the expression for this matrix entry
      ex = "0.0"
      for sym = 0:max_symbol_idx(l1, l2)
         # expression for the new entry
         # ex1 = Gsym(l1, l2, m1, m2, sym)
         ex1 = tbl[_lookupkey(l1,l2,m1,m2,sym)]
         # extression for the bond integral V
         V_idx = bondintegral_index(l1, l2, sym)
         ex = "$ex + ($ex1) * V[$V_idx]"
      end
      ex_assign = " g[$I1, $I2] = $ex "
      push!(code, Calculus.simplify(Meta.parse(ex_assign)))

      # symmetric part (if not on the diagonal)
      if I2 != I1
         # sign modification
         sig = sksign(l2, l1) * signmod(SGN)
         push!(code, :( g[$I2, $I1] = $sig * g[$I1, $I2] ))
      end
   end
   quote
      $(Expr(:block, code...))
      return g
   end
end


# ----------------------------------------------------------------------
#   NEW IMPLEMENTATION

@generated function sk_gen(::SKBond{O1, O2, SYM}, φ, θ) where {O1, O2, SYM}
   # get the SK expressions table
   l1, l2 = get_l(Val{O1}()), get_l(Val{O2}())
   tbl = sk_table(max(l1,l2))
   expr_str = "SMatrix{$(2*l1+1),$(2*l2+1)}("
   for m2 = -l2:l2, m1 = -l1:l1
      ex = tbl[_lookupkey(l1, l2, m1, m2, get_bidx(Val{SYM}()))]
      expr_str *= ex * ", "
   end
   code = Meta.parse(expr_str[1:end-2] * ")")
   quote
      $code
   end
end

end
