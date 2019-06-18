
@testset "Slater-Koster-Core" begin

@info("Slater Koster Core Tests...")
using SlaterKoster, Test, LinearAlgebra
import SlaterKoster.CodeGeneration
using SlaterKoster: SKH, sk2cart, cart2sk, allbonds, nbonds

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


@info("New implementation: sp")
orbitals = [sko"s", sko"p"]
bonds = [skb"ssσ", skb"spσ", skb"ppσ", skb"ppπ"]
println(@test bonds == allbonds(orbitals))
# H = SKH(orbitals, bonds)   # TODO: need to fix this test somehow!!! (or drop it)
# println(@test H == SKH("sp"))
H = skh"sp"
for n = 1:5
   V = rand(length(bonds))
   U = rand(3) .- 0.5
   U /= norm(U)
   Hnew = sk2cart(H, U, V)
   Hold = SK.OldSK.sk4!(U, V, zeros(4,4))
   perm = [1,3,4,2]
   println(@test Hnew ≈ Hold[perm, perm])
end

@info("New implementation: spd")
orbitals = [sko"s", sko"p", sko"d"]
H = SKH(orbitals) # generate bonds automatically
println(@test H == SKH("spd"))
H = SKH("spd")
for n = 1:5
   V = rand(nbonds(H))
   U = rand(3) .- 0.5
   U /= norm(U)
   Hnew = sk2cart(H, U, V)
   Hold = SK.OldSK.sk9!(U, V, zeros(9,9))
   perm = [1,3,4,2,5,6,9,7,8]
   println(@test Hnew ≈ Hold[perm, perm])
end

@info("sk2cart to cart2sk : sp")
orbitals = [sko"s", sko"p"]
H = SKH(orbitals) # generate bonds automatically
println(@test bonds == allbonds(orbitals))
H = SKH("sp")
for n = 1:5
   V = rand(nbonds(H))
   U = rand(3) .- 0.5
   U /= norm(U)
   Hnew = sk2cart(H, U, V)
   Vnew = cart2sk(H, U, Hnew)
   println(@test V ≈ Vnew)
end

@info("sk2cart to cart2sk : spd")
orbitals = [sko"s", sko"p", sko"d"]
H = SKH(orbitals) # generate bonds automatically
println(@test H == SKH("spd"))
H = SKH("spd")
for n = 1:5
   V = rand(nbonds(H))
   U = rand(3) .- 0.5
   U /= norm(U)
   Hnew = sk2cart(H, U, V)
   Vnew = cart2sk(H, U, Hnew)
   println(@test V ≈ Vnew)
end

end
