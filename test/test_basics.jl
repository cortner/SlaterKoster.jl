
@testset "SlaterKoster Basics" begin

   using SlaterKoster: Orbital, parse_orbital

   @test parse_orbital("s") == Orbital(:s)


end


using  SlaterKoster, Test
using SlaterKoster: Orbital, parse_orbital, get_sym, Bond, parse_bond

function test_fail(f)
   try
      f()
      @info("We should have thrown an error here")
      return false
   catch
      return true
   end
end

@test parse_orbital("s") == parse_orbital(:s)
@test parse_orbital("2p") != parse_orbital(:_2p)
@test test_fail( _ -> parse_orbital("d") == parse_orbital("3d"))
@test test_fail( _ -> parse_orbital("x") )
o = Orbital("3d")

b = parse_bond("2s2pσ")
@test
