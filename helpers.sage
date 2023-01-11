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
    coef = G.order() / (G.centralizer(a).order() * G.centralizer(b).order())
    return coef * sum


def simple_objects(G):
    reps = G.conjugacy_classes_representatives()
    simples = []
    for g in reps:
        irr_chs = G.centralizer(g).irreducible_characters()
        g_simples = [(g, irr_ch) for irr_ch in irr_chs]
        simples += g_simples
    return simples


def s_matrix(G):
    simples = simple_objects(G)
    rows = []
    for simple in simples:
        row = [s_matrix_entry(G, simple, sp) for sp in simples]
        rows.append(row)
    return Matrix(rows)
