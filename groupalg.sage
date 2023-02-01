"""
Minimal implementation of the group algebra
"""

from sage.categories.all import Algebras, AlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.rings.integer_ring import ZZ
from sage.misc.misc import inject_variable

class GAlg(CombinatorialFreeModule):
    def __init__(self, G, prefix="s"):
        self._G = G
        self._prefix = prefix
        self._names = {G.one() : "%s0"%self._prefix}
        cat = AlgebrasWithBasis(ZZ).Subobjects()
        CombinatorialFreeModule.__init__(self, ZZ, self._G, category=cat)
        count = 1
        for g in G:
            if g == G.one():
                continue
            name = "%s%s" % (prefix, count)
            self._names[g] = name
            count += 1
        for g in self._G:
            inject_variable(self._names[g],self.monomial(g))

    def _repr_(self):
        return "Group algebra of %s"%self._G

    def __call__(self, *args):
        if len(args) > 1:
            args = (args,)
        return super(GAlg, self).__call__(*args)

    def _element_constructor(self, g):
        return self.monomial(g)

    def product_on_basis(self, a, b):
        d = {a*b: 1}
        return self._from_dict(d)

    def _repr_term(self, t):
        return self._names[t]

    class Element(CombinatorialFreeModule.Element):
        def g(self):
           """
           Return the underlying group element of a monomial element
           """
           return self.support_of_term()

        def order(self):
            """
            Return the order of a monomial element
            """
            return self.support_of_term().order()
