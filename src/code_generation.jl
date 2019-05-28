
__precompile__(false)
module CodeGeneration

using PyCall, Calculus

using SlaterKoster: max_symbol, bondintegral_index, orbital_index


# load the python part of the code generation
const _pysys = pyimport("sys")
push!(_pysys."path", @__DIR__())
const _codegen = pyimport("codegen")


function signmatrix(l_1, l_2)
   S = fill(1, (2*l_1+1, 2*l_2+1))
   for i = 1:2*l_1+1, j = 1:2*l_2+1
      if i - l_1 > 0 # && (i - l_1) % 2 == 1
         S[i,j] *= -1
      end
      if j - l_2 > 0 # && (j - l_2) % 2 == 1
         S[i,j] *= -1
      end
   end
   return S
end



"""
## INPUTS

'l_1': Angular quantum number of atom 1.
'l_2': Angular quantum number of atom 2.
     l        : 0  1  2  3  4
     orbitals : s  p  d  f  g  h  j  k ...
'm_1': Magnetic quantum number of atom 1.
'm_2': Magnetic quantum number of atom 2.
     m : sub-orbitals for the ls

'M': A float which represents the type of symmetry bond we are considering,
        M | 0 1 2 3 4
    symbol| σ π δ φ γ

## RETURNS

This function returns the relevant coefficient that should be are used in writing the
Slater-Koster transformations.
"""
function Gsym(l_1, l_2, m_1, m_2, M::Integer)
   str = _codegen.Gsym(l_1, l_2, m_1, m_2, M)
   # replace some symbols to turn python into julia code
   str = replace(str, "**" => "^")
   # α : Azithmuthal coordinate (normally φ)
   str = replace(str, "alpha" => "φ")   # the python code uses alpha here
   # β : Polar angle (normally θ)
   str = replace(str, "beta" => "θ")   # the python code uses alpha here
   return str
end

const sigs = [1 1 1 1; -1 1 1 1; -1 1 1 1; -1 1 1 1]

@generated function sk_gen!(g, ::Val{L}, V, φ, θ) where {L}
   code = Expr[]
   for l1 = 0:L, l2 = 0:L, m1 = -l1:l1, m2 = -l2:l2
      I1 = orbital_index(l1, m1)
      I2 = orbital_index(l2, m2)
      sig = sigs[I1, I2]
      ex = "0.0"
      for sym = 0:max_symbol(l1, l2)
         # expression for the new entry
         ex1 = Gsym(l1, l2, m1, m2, sym)
         # extression for the bond integral V
         V_idx = bondintegral_index(l1, l2, sym)
         ex = "$ex + ($ex1) * V[$V_idx]"
      end
      ex_assign = " g[$I1, $I2] = $sig * ($ex) "
      push!(code, Calculus.simplify(Meta.parse(ex_assign)))
   end
   quote
      $(Expr(:block, code...))
      return g
   end
end


end
