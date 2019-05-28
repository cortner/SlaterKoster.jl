
# part of SlaterKoster.jl

max_symbol(l1,l2) = min(l1, l2)

MSym = Dict( 0 => "σ", 1 => "π", 2 => "δ", 3 => "φ", 5 => "γ" )
LSym = Dict( 0 => "s", 1 => "p", 2 => "d", 3 => "f", 4 => "g" )

"""
`bondintegral_index` :

 1 ssσ  -> (0, 0, 0) = (l1, l2, M) -> 1
 2 psσ
 3 ppσ
 4 ppπ
 5 dsσ
 6 dpσ
 7 dpπ
 8 ddσ
 9 ddπ
10 ddδ

TODO: there must be a cleverer implementation of this!
"""
function bondintegral_index(l1, l2, M)
   if l1 < l2
      l1, l2 = l2, l1
   end
   @assert M <= max_symbol(l1, l2)

   idx = 0
   for _l1 = 0:l1, _l2 = 0:_l1, _M = 0:max_symbol(_l1, _l2)
      idx += 1
      if (_l1, _l2, _M) == (l1, l2, M)
         return idx
      end
   end
   @error("we shouldn't be here!")
end

"""
     l  m
s    0  0    -> 1
pz   1  -1   -> 2
py   1  0    -> 3
px   1   1   -> 4
...

TODO: there must be a cleverer implementation of this!
"""
function orbital_index(l, m)
   @assert abs(m) <= l
   idx = 0
   for _l = 0:l, _m = -_l:_l
      idx += 1
      if (_l, _m) == (l, m)
         return idx
      end
   end
   @error("we shouldn't be here!")
end
