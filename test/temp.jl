
using SlaterKoster, Test, LinearAlgebra, Calculus
import SlaterKoster.CodeGeneration

SK = SlaterKoster
CG = SlaterKoster.CodeGeneration

@show CG.Gsym(1, 0, 1, 0, 0)
CG.Gsym(0,0,0,0,0)
CG.Gsym(0,1,0,-1,0)

@show CG.Gsym(1,1,1,1,0)

function expr(l1, l2, m1, m2)
   ex = :(0.0)
   for M = 0:min(l1,l2)
      ex1 = Meta.parse(CG.Gsym(l1, l2, m1, m2, M))
      # Vsym = Meta.parse("V$(SK.LSym[l1])$(SK.LSym[l2])$(SK.MSym[M])")
      Vsym = Meta.parse("V[$(SK.bondintegral_index(l1, l2, M))]")
      ex = Calculus.simplify(:($(ex)  + $(ex1) * $(Vsym)))
   end
   I1 = SK.orbital_index(l1, m1)
   I2 = SK.orbital_index(l2, m2)
   println("mat[$I1, $I2]")
   println("   ", ex)
end

for l1 = 0:1, l2 = l1:1 # s,p,d
   for m1=-l1:l1, m2=-l2:l2
      ex = :(0.0)
      for M = 0:min(l1,l2)
         ex1 = Meta.parse(CG.Gsym(l1, l2, m1, m2, M))
         # Vsym = Meta.parse("V$(SK.LSym[l1])$(SK.LSym[l2])$(SK.MSym[M])")
         Vsym = Meta.parse("V[$(SK.bondintegral_index(l1, l2, M))]")
         ex = Calculus.simplify(:($(ex)  + $(ex1) * $(Vsym)))
      end
      I1 = SK.orbital_index(l1, m1)
      I2 = SK.orbital_index(l2, m2)
      println("mat[$I1, $I2]")
      println("   ", ex)
   end
end


# joinpath(@__DIR__(), "..", "src")
# using PyCall
# pysys = pyimport("sys")
# push!(pysys."path", joinpath(@__DIR__(), "..", "src"))
# codegen = pyimport("codegen")
# @show codegen.Gsym(1,1,0,0,0)
# # ex = Meta.parse(replace(codegen.Gsym(1,1,0,0,0), "**" => "^"))
# # ex1 = :( $ex + $ex )
# # f = eval(:(beta -> $ex1))
# # f(0.0)
