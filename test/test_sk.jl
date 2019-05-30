
@testset "Slater-Koster-Core" begin

@info("Slater Koster Core Tests...")
using SlaterKoster, Test, LinearAlgebra
import SlaterKoster.CodeGeneration

SK = SlaterKoster
CG = SlaterKoster.CodeGeneration

@info("carttospher test")
for n = 1:10
   U = rand(3) .- 0.5
   U = U / norm(U)
   φ, θ = SK.carttospher(U...)
   r, α, β = CG._codegen.carttospher(U...)
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
   Hold = zeros(9,9)
   Hnew = zeros(9,9)
   SK.OldSK.sk9!(U, V, Hold)
   SK._sk!(Hnew, Val(2), U, V)
   perm = [1,3,4,2,5,6,9,7,8]
   print((@test Hold[perm, perm] ≈ Hnew), " ")
   # @test sort(eigvals(Hold)) ≈ sort(eigvals(Hnew))
end
println()

end
