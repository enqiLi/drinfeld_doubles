#!/usr/bin/env sage


from helpers import Double
from sage.all import *
import sys
import os
from timeit import default_timer as timer

# attach('./helpers.sage')


cyclics = [CyclicPermutationGroup(n) for n in range(2, 7)]
dihedrals = [DihedralGroup(n) for n in range(4, 30)]
alternatings = [AlternatingGroup(n) for n in range(4, 6)]
symmetrics = [SymmetricGroup(n) for n in range(3, 6)]
groups = cyclics + dihedrals + alternatings + symmetrics + [KleinFourGroup()]

for G in groups:
    start = timer()
    D = Double(G)
    B = D.basis()
    Ns = [D.N_ijk(x, y, z) for x in B
          for y in B for z in B]
    if all(c <= 1 for c in Ns):
        print(G, " is multiplicity free.")
    else:
        print(G, " is NOT multiplicity free.")
    end = timer()
    print("%s seconds is used for the above group. It has %s simple objects.\n" % (end - start, len(B)))
