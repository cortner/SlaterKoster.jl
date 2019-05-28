using SlaterKoster, Test, LinearAlgebra, Calculus
import SlaterKoster.CodeGeneration

SK = SlaterKoster
CG = SlaterKoster.CodeGeneration


for n = 1:10
   U = rand(3) .- 0.5
   U = U / norm(U)
   φ, θ = SK.carttospher(U...)
   V = [sin(φ)*sin(θ), cos(φ)*sin(θ), cos(θ)]
   @test U ≈ V
end

for n = 1:10
   V = rand(4) .- 0.5
   U = rand(3) .- 0.5
   U = U / norm(U)
   Hold = zeros(4,4)
   Hnew = zeros(4,4)
   SK.OldSK.sk4!(U, V, Hold)
   SK._sk!(Hnew, Val(1), U, V)
   perm = [1,2,4,3]
   @test Hold[perm, perm] ≈ Hnew
   # @test sort(eigvals(Hold)) ≈ sort(eigvals(Hnew))
end



# @testset "SlaterKoster" begin
#
#    include("test_kwon.jl")
#
# end
