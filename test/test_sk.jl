
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
   Hold = SK.OldSK.sk9!(U, V, zeros(9,9))
   Hnew = SK.sk(Hnew, Val(2), U, V)
   perm = [1,3,4,2,5,6,9,7,8]
   print((@test Hold[perm, perm] ≈ Hnew), " ")
   # @test sort(eigvals(Hold)) ≈ sort(eigvals(Hnew))
end
println()

@info("test the FHI-aims format")
set1 = JSON.parsefile(@__DIR__() * "/data/sp_o3_offsite_data.json"
set2 = JSON.parsefile(@__DIR__() * "/data/spdf_au2_offsite_data.json"
for set in (set1, set2)

end

end
