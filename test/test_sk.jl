
@testset "Slater-Koster-Core" begin

@info("Slater Koster Core Tests...")
using SlaterKoster, Test, LinearAlgebra
import SlaterKoster.CodeGeneration
using SlaterKoster: SKH, sk, _sk!, allbonds, nbonds

SK = SlaterKoster
CG = SlaterKoster.CodeGeneration

# temporarily keep this to test against python code
using PyCall
_pysys = pyimport("sys")
push!(_pysys."path", @__DIR__()[1:end-4]*"src")
_codegen = pyimport("codegen")

@info("carttospher test")
for n = 1:5
   U = rand(3) .- 0.5
   U = U / norm(U)
   φ, θ = SK.carttospher(U...)
   r, α, β = _codegen.carttospher(U...)
   print((@test φ ≈ α), " ")
   print((@test θ ≈ β), " ")
   V = [cos(φ)*sin(θ), sin(φ)*sin(θ), cos(θ)]
   print((@test U ≈ V), " ")
end
println()

@info("sk! test (the first sk! call is slow!)")
for n = 1:10
   V = rand(10) .- 0.5
   U = rand(3) .- 0.5
   U = U / norm(U)
   Hold = SK.OldSK.sk9!(U, V, zeros(9,9))
   Hnew = SK._sk!(zeros(9,9), Val(2), U, V)
   perm = [1,3,4,2,5,6,9,7,8]
   print((@test Hold[perm, perm] ≈ Hnew), " ")
end
println()


@info("New implementation: sp")
orbitals = [sko"s", sko"p"]
bonds = [skb"ssσ", skb"spσ", skb"ppσ", skb"ppπ"]
println(@test bonds == allbonds(orbitals))
H = SKH(orbitals, bonds)
println(@test H == SKH("sp"))
for n = 1:5
   V = rand(length(bonds))
   U = rand(3) .- 0.5
   U /= norm(U)
   Hnew = sk(H, U, V)
   Hold = SK.OldSK.sk4!(U, V, zeros(4,4))
   perm = [1,3,4,2]
   println(@test Hnew ≈ Hold[perm, perm])
end

@info("New implementation: spd")
orbitals = [sko"s", sko"p", sko"d"]
H = SKH(orbitals) # generate bonds automatically
println(@test H == SKH("spd"))
for n = 1:5
   V = rand(nbonds(H))
   U = rand(3) .- 0.5
   U /= norm(U)
   Hnew = sk(H, U, V)
   Hold = SK.OldSK.sk9!(U, V, zeros(9,9))
   perm = [1,3,4,2,5,6,9,7,8]
   println(@test Hnew ≈ Hold[perm, perm])
end

end
