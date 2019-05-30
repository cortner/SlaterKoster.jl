using SlaterKoster, JuLIP, Test, LinearAlgebra

@testset "Kwon Test" begin
   @info("Kwon Tests...")
   kwon = KwonHamiltonian()
   at = bulk(:Si, cubic = true, pbc = true) * 4
   H, M = hamiltonian(kwon, at)
println(@test M == Matrix((1.0+0.0im)*I, size(M)))

# -----------------------------------------------------------------------------
# visual check that the spectrum looks ok - need to make this a proper
# test set
# -----------------------------------------------------------------------------
# E = eigvals(H)
# using Plots
# plot(E, lw=0, ms=5, label = "")
# eF = 0.5 * (E[length(E)÷2] + E[length(E)÷2+1])
# hline!([eF], style=:dot, c=:black, label ="")
# -----------------------------------------------------------------------------
end
