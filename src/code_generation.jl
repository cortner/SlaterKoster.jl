
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

# TODO: switch this to SymPy.jl ???

sympy = pyimport("sympy")
sp_spin = pyimport("sympy.physics.quantum.spin")

const spsin = sympy.sin
const spcos = sympy.cos
const Rotation = sp_spin.Rotation

_lookupkey(l1, l2, m1, m2, sym) = "$l1,$l2,$m1,$m2,$sym"
_fname_sktable() = joinpath(@__DIR__(), "sktable.json")

"""
`small_d(l, m, M, θ)` : Numerical calculation for Wigner small d-function

"""
function small_d(l, m, mp, θ)
    fc1 = factorial(l+m)
    fc2 = factorial(l-m)
    fc3 = factorial(l+mp)
    fc4 = factorial(l-mp)
    fcm1 = sqrt(fc1 * fc2 * fc3 * fc4)
  
    cosb = cos(θ / 2.0)
    sinb = sin(θ / 2.0)

    p = m - mp
    lo = max(0,p)
    hi = min(0,l+m,l-mp)
   
    rtn_sum = 0.0
    for k = lo:hi
       fc5 = factorial(k)
       fc6 = factorial(l+m-k)
       fc7 = factorial(l-mp-k)
       fc8 = factorial(k-p)
       fcm2 = fc5 * fc6 * fc7 * fc8
       pow1 = 2 * l - 2 * k + p
       pow2 = 2 * k - p
       rtn_sum += (-1)^(k+p) * cosb^pow1 * sinb^pow2 / fcm2
    end
    rtn_sum *= fcm1

    return rtn_sum
end


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
   tbl = Dict{String, Any}("L" => L)
   for l1 = 0:L, l2=0:L, m1=-l1:l1, m2=-l2:l2
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

   #expr = expr.simplify()

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

tau(m) = (m >= 0) ? 1 : 0

AA(m, al) = ( (m == 0) ? sqrt(2)/2 :
            (-1)^m * (tau(m) * cos(abs(m)*al) + tau(-m) * sin(abs(m)*al)) )

BB(m, al) = ( (m == 0) ? 0 :
            (-1)^m * (tau(-m) * cos(abs(m)*al) - tau(m) * sin(abs(m)*al) ) )

SS(l, m, sym, al, be) = AA(m, al) * ( (-1)^sym * small_d(l, abs(m),  sym, be)
                                               + small_d(l, abs(m), -sym, be) )

TT(l, m, sym, al, be) = BB(m, al) * ( (-1)^sym * small_d(l, abs(m),  sym, be)
                                               - small_d(l, abs(m), -sym, be) )

function Gnum(l1, l2, m1, m2, sym, al, be)
  
   rtn = 0.0

   if sym == 0
      rtn += (2 * AA(m1, al) * small_d(l1, abs(m1), 0, be)
                * AA(m2, al) * small_d(l2, abs(m2), 0, be) )
   else
      rtn += (   SS(l1, m1, sym, al, be) * SS(l2, m2, sym, al, be)
               + TT(l1, m1, sym, al, be) * TT(l2, m2, sym, al, be) )
   end

   return rtn
end

function sk_num(::SKBond{O1, O2, SYM}, φ, θ) where {O1, O2, SYM}
   l1, l2 = get_l(Val{O1}()), get_l(Val{O2}())
   E = zeros(2*l1+1,2*l2+1)
   b_l = get_bidx(Val{SYM}())
   for (mj, m2) = enumerate(-l2:l2), (mi, m1) = enumerate(-l1:l1)
      E[mi, mj] = Gnum(l1, l2, m1, m2, b_l, φ, θ) 
   end
   return E
end

end

# TODO: create in-place type-stable versions of sk_gen
#       taking l, m as input 
# _parse_key(k) = eval(Meta.parse("[" * k * "]"))
#
# # generate code for sk blocks
# # require at least spd
# tbl = sk_table(2)
# lmax = 2
# for (key, expr) in tbl
#    lmax = max(lmax, _parse_key(key)[1])
# end
# _sk_f_tbl =
