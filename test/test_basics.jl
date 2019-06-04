
@testset "SlaterKoster Basics" begin

   using  SlaterKoster, Test
   using SlaterKoster:  get_l, get_bidx

   function test_fail(f)
      try
         f()
         @info("We should have thrown an error here")
         return false
      catch
         return true
      end
   end

   @info("Testing `get_l` for orbitals")
   println(@test get_l(Val{:s}()) == 0)
   println(@test get_l(Val{:p}()) == 1)
   println(@test get_l(Val{:d}()) == 2)
   println(@test get_l(Val{:f}()) == 3)
   println(@test get_l(Val{:g}()) == 4)

   @info("Testing `get_l, get_bidx` for bonds")
   println(@test get_l(skb"spσ") == (0, 1))
   println(@test get_bidx(skb"spσ") == 0)
   println(@test get_l(skb"2p2pπ") == (1, 1))
   println(@test get_bidx(skb"2p2pπ") == 1)
   println(@test get_l(skb"1s2pσ") == (0,1))
   println(@test get_bidx(skb"1s2pσ") == 0)
end
