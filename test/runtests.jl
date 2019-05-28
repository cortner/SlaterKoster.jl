
using SlaterKoster, Test
import SlaterKoster.CodeGeneration

SK = SlaterKoster
CG = SlaterKoster.CodeGeneration

@show CG.Gsym(1, 1, 1, 1, 1)
CG.Gsym(0,0,0,0,0)
CG.Gsym(0,1,0,-1,0)

for l1 = 0:1, l2 = l1:1 # s,p,d
   for m1=-l1:l1, m2=-l2:l2
      ex = :(0.0)
      for M = 0:min(l1,l2)
         ex1 = CG.Gsym(l1, l2, m1, m2, M)
         Vsym = Meta.parse("V$(SK.LSym[l1])$(SK.LSym[l2])$(SK.MSym[M])")
         ex = Calculus.simplify(:($(ex)  + $(ex1) * $(Vsym)))
      end
      # println("$(SK.LSym[l1]), $l2, $m1, $m2")
      println("   ", ex)
   end
end





# @testset "SlaterKoster" begin
#
#    include("test_kwon.jl")
#
# end



# joinpath(@__DIR__(), "..", "src")
# using PyCall
# pysys = pyimport("sys")
# push!(pysys."path", joinpath(@__DIR__(), "..", "src"))
# codegen = pyimport("codegen")
# @show codegen.Gsym(1,1,0,0,0)
# ex = Meta.parse(replace(codegen.Gsym(1,1,0,0,0), "**" => "^"))
# ex1 = :( $ex + $ex )
# f = eval(:(beta -> $ex1))
# f(0.0)
