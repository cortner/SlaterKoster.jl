

@testset "Kwon Test" begin

using SlaterKoster, JuLIP, Test, LinearAlgebra
OldSK = SlaterKoster.OldSK

using SlaterKoster: assembleRI

@info("Old Kwon Tests... just checking that it still runs in principle")
oldkwon = OldSK.KwonHamiltonian()
at = bulk(:Si, cubic = true, pbc = false) * 2
H, M = OldSK.hamiltonian(oldkwon, at)
println(@test M == Matrix((1.0+0.0im)*I, size(M)))

@info("Test Old Kwon versus new Kwon hamiltonians")
ats = [ rattle!(bulk(:Si, cubic = true, pbc = false) * 2, 0.1),
        rattle!(bulk(:Si, cubic = true, pbc = true) * 2, 0.1),
        rattle!(bulk(:Si, cubic = false, pbc = true) * 3, 0.1) ]
kwon = KwonTB()

for at in ats
   Hold, Mold = OldSK.hamiltonian(oldkwon, at)
   # for i = 1:size(Hold,1); Hold[i,i] = 0; end
   eigold = eigvals(Hold)

   Hri = assembleRI(at, kwon)
   Hc = Matrix(Hri)
   eignew = eigvals(Symmetric(Hc))

   println(@test sort(eigold) ≈ sort(eignew))
end

end


# # -----------------------------------------------------------------------------
# # visual check that the spectrum looks ok - need to make this a proper
# # test set
# # -----------------------------------------------------------------------------
# # E = eigvals(H)
# # using Plots
# # plot(E, lw=0, ms=5, label = "")
# # eF = 0.5 * (E[length(E)÷2] + E[length(E)÷2+1])
# # hline!([eF], style=:dot, c=:black, label ="")
# # -----------------------------------------------------------------------------
