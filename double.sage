from sage.categories.all import Algebras, AlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.misc.misc import inject_variable
from sage.rings.integer_ring import ZZ

def summands(G, a, b):
    return [g for g in G.list() if a*g*b*g ^ (-1) == g*b*g ^ (-1)*a]


def s_matrix_entry(G, simple_1, simple_2):
    '''
    simple_1 and simple_2 are simple objects, which are tuples
    (g, chi) where chi is an irreducible character of the centralizer
    of g in G.
    '''
    sum = 0
    a, chi_1 = simple_1
    b, chi_2 = simple_2
    for g in summands(G, a, b):
        sum += chi_1(g*b*g ^ (-1)) * chi_2(g ^ (-1)*a*g)
    coef = 1 / (G.centralizer(a).order() * G.centralizer(b).order())
    return coef * sum


def simple_objects(G):
    reps = G.conjugacy_classes_representatives()
    simples = []
    for g in reps:
        irr_chs = G.centralizer(g).irreducible_characters()
        g_simples = [(g, irr_ch) for irr_ch in irr_chs]
        simples += g_simples
    return simples


class Simple(CombinatorialFreeModule.Element):
    def __init__(self, group, ccrep, irr, name, parent):
        """
        A class for simple objects in a Drinfeld double

        Input group - a finite group
        ccrep - a conjugacy class representative
        irr - an irrep of the centralizer of g
        """
        self._group = group
        self._cc = ccrep
        self._irr = irr
        self._name = name
        self._parent = parent

    def __repr__(self):
        return self._name

    def irr(self):
        return self._irr

    def g(self):
        return self._cc

    def group(self):
        return self._group

    def parent(self):
        return self._parent

    def dual(self):
        return self._parent._conjugates[self]

    def twist(self):
        return self._parent.twist(self)

def s_matrix(G):

    simples = simple_objects(G)
    rows = []
    for simple in simples:
        row = [s_matrix_entry(G, simple, sp) for sp in simples]
        rows.append(row)
    return Matrix(rows)

class Double(CombinatorialFreeModule):
    """
    A class for the Drinfeld double of a finite group
    """
    def __init__(self, group, prefix="s"):
        self._group = group
        self._prefix = prefix
        self._simples = []
        count = 0
        
        for g in group.conjugacy_classes_representatives():
            for irr in group.centralizer(g).irreducible_characters():
                name = "%s%s" % (prefix, count)
                obj = Simple(group, g, irr, name, self)
                self._simples.append(obj)
                inject_variable(name, obj)
                count += 1
        # This doesn't seem to work
        cat = AlgebrasWithBasis(ZZ).Subobjects()
        CombinatorialFreeModule.__init__(self, ZZ, basis_keys = self._simples, element_class=Simple, category=cat)
        self._smatrix = {(i, j): self._s(i, j)
                         for i in self._simples for j in self._simples}
        self._conjugates = {i: j for i in self._simples for j in self._simples if self.N_ijk(i,j,self._simples[0]) == 1}

    def __repr__(self):
        return "Double of %s" % self._group
    
    def __call__(self, *args):
        # Newly added. Not sure about this.
        if len(args) > 1:
            args = (args,)
        return super(Double, self).__call__(*args)
    
    def _element_constructor(self, g, irr, name):
        # Newly added
        return Simple(self._group, g, irr, name, self)

    def _s(self, i, j):
        """
        The s-matrix entry corresponding to simple objects a, b
        """
        a, chi_1 = i.g(), i.irr()
        b, chi_2 = j.g(), j.irr()
        G = self._group
        ret = sum(chi_1(g*b*g ^ (-1)) * chi_2(g ^ (-1)*a*g)
                  for g in summands(G, a, b))
        coef = 1 / (G.centralizer(a).order() * G.centralizer(b).order())
        return coef * ret

    def _repr_term(self, t):
        # Not the right way, but creating a name list causes problems
        return t._name
    
    def product_on_basis(self, a, b):
        # Newly added
        d = {c: self.N_ijk(a, b, c.dual()) for c in self.basis()}
        return self._from_dict(d)

    def s(self, i, j):
        return self._smatrix[(i, j)]

    def s0(self):
        return self._simples[0]

    def basis(self):
        return self._simples

    def N_ijk(self, i, j, k):
        sz = self.s0()
        return ZZ(sum(self.s(i, r) * self.s(j, r) * self.s(k, r)/self.s(sz, r) for r in self.basis()))

    def twist(self, i):
        """
        The twist of the simple object i
        """
        return i.irr()(i.g()) / i.irr()(self._group("()"))

    def p_plus(self):
        return self._group.order()

    def p_minus(self):
        return self._group.order()
