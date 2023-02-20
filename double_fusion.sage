"""
The Fusion Ring of the Drinfeld Double of a Finite Group
"""
# ****************************************************************************
#  Copyright (C) 2023 Wenqi Li
#                     Daniel Bump <bump at match.stanford.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.categories.all import Algebras, AlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.rings.integer_ring import ZZ
from sage.misc.misc import inject_variable
from sage.misc.cachefunc import cached_method

class FusionDouble(CombinatorialFreeModule):
    r"""
    This constructs the Fusion Ring of the modular
    tensor category of modules over the Drinfeld
    Double of a finite group. Usage may be similar
    Sage FusionRings.

    INPUT:

    - ``G`` -- a finite group
    - ``prefix`` (optional: defaults to 's') a prefix
      for the names of simple objects.

    EXAMPLES ::
        sage: G = DihedralGroup(5)
        sage: H = FusionDouble(G)
        sage: H.basis()
        Finite family {0: s0, 1: s1, 2: s2, 3: s3, 4: s4, 5: s5, 6: s6, 7: s7, 8: s8, 9: s9, 10: s10, 11: s11, 12: s12, 13: s13, 14: s14, 15: s15}
        sage: for x in H.basis():
        ....:     print ("%s : %s"%(x,x^2))
        ....: 
        s0 : s0
        s1 : s0
        s2 : s0 + s1 + s3
        s3 : s0 + s1 + s2
        s4 : s0 + s2 + s3 + s6 + s7 + s8 + s9 + s10 + s11 + s12 + s13 + s14 + s15
        s5 : s0 + s2 + s3 + s6 + s7 + s8 + s9 + s10 + s11 + s12 + s13 + s14 + s15
        s6 : s0 + s1 + s11
        s7 : s0 + s1 + s13
        s8 : s0 + s1 + s15
        s9 : s0 + s1 + s12
        s10 : s0 + s1 + s14
        s11 : s0 + s1 + s6
        s12 : s0 + s1 + s9
        s13 : s0 + s1 + s7
        s14 : s0 + s1 + s10
        s15 : s0 + s1 + s8
        sage: s4*s5
        s1 + s2 + s3 + s6 + s7 + s8 + s9 + s10 + s11 + s12 + s13 + s14 + s15
        sage: s8.twist()
        zeta5^3
    """
    def __init__(self, G, prefix="s"):
        self._G = G
        self._prefix = prefix
        self._names = {}
        self._elt = {}
        self._chi = {}
        count = 0
        for g in G.conjugacy_classes_representatives():
            for chi in G.centralizer(g).irreducible_characters():
                self._names[count] = "%s%s"%(prefix, count)
                self._elt[count] = g
                self._chi[count] = chi
                count += 1
        self._rank = count
        cat = AlgebrasWithBasis(ZZ).Subobjects()
        CombinatorialFreeModule.__init__(self, ZZ, [k for k in self._names], category=cat)
        for i in range(self._rank):
            inject_variable(self._names[i],self.monomial(i))

    def _repr_(self):
        return "The Fusion Ring of the Drinfeld Double of %s"%self._G

    def __call__(self, *args):
        if len(args) > 1:
            args = (args,)
        return super(GAlg, self).__call__(*args)

    def _element_constructor(self, k):
        return self.monomial(k)

    @cached_method
    def s_ij(self, i, j):
        sum = 0
        [i, j] = [x.support_of_term() for x in [i,j]]
        [a, chi_1] = [self._elt[i], self._chi[i]]
        [b, chi_2] = [self._elt[j], self._chi[j]]
        G = self._G
        for g in G:
            if a*g*b*g.inverse() == g*b*g.inverse()*a:
                sum += chi_1(g*b*g ^ (-1)) * chi_2(g ^ (-1)*a*g)
        coef = 1 / (G.centralizer(a).order() * G.centralizer(b).order())
        return coef * sum

    @cached_method
    def N_ijk(self, i, j, k):
        """
        The symmetric invariant of three simple objects,
        this returns the dimension of

        .. MATH::
           Hom(i \\otimes j\\otimes k, s0)

        where `s_0` is the unit element (assuming prefix='s').
        """
        sz = self.one()
        return ZZ(sum(self.s_ij(i, r) * self.s_ij(j, r) * self.s_ij(k, r)/self.s_ij(sz, r) for r in self.basis()))
    

    def Nk_ij(self, i, j, k):
        r"""
        Returns the fusion coefficient `N^k_{ij}`
        """
        return self.N_ijk(i, j, self.dual(k))

    def one(self):
        return self.basis()[0]

    def dual(self,i):
        sz = self.one()
        for j in self.basis():
            if self.N_ijk(i,j,sz) > 0:
                return j

    def product_on_basis(self, a, b):
        d = {k.support_of_term() : self.N_ijk(self.monomial(a),self.monomial(b),self.dual(k)) for k in self.basis()}
        return self._from_dict(d)

    def _repr_term(self, t):
        return self._names[t]

    class Element(CombinatorialFreeModule.Element):
        def g(self):
           return self.parent()._elt[self.support_of_term()]

        def chi(self):
           return self.parent()._chi[self.support_of_term()]
    
        def twist(self):
            """
            The twist of the simple object i
            """
            return self.chi()(self.g()) / self.chi()(self.parent()._G.one())
