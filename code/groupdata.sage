"""
Computes basic data associated to group Gamma0(N)

The key tool is group_data(N), which returns a dictionary with the following
keys and values:

    ### group_data(N) dictionary ###

    KEY        : VALUE DESCRIPTION
    -----------:------------------
    coset_reps : List of coset representatives. Each element in the list is a sage matrix.
    cusps      : List of cusps. Infinity is the first cusp.
    pvertices  : Set of parabolic vertices.
    level      : Integer level
    Udata      : Udict (see below)
    Cdata      : Cdict (see below)

The two subdictionaries Udata and Cdata store data associated to the parabolic
vertices and cusps. The format of these dictionaries is the same except for the
keys.

    ### Udict dictionary ###

    KEYS are parabolic vertices in fundamental domain.
    Values are tuples (pvertex, cusp, U_ell, sigma_ell), where
        pvertex is the parabolic vertex,
        cusp is the equivalent cusp class for `pvertex`,
        U_ell is a Gamma0(N) matrix taking `pvertex` to `cusp`, and
        sigma_ell is the cusp normalizing map for `cusp`.

    Note: the first entry is the same as the key.


    ### Cdict dictionary ###

    KEYS are cusps
    Values are tuples (cusp, cusp, U_ell, sigma_ell), where
        cusp is the cusp,
        the second cusp is also the cusp,
        U_ell is the identity matrix, and
        sigma_ell is the cusp normalizing map for `cusp`.

    Note: the structure of Cdict is to mirror Udict, even though much of the
    data is redundant.

"""

def check_level(N):
    """
    Checks that N = 2^r p_1 ... p_n, r <= 3 and p_j distinct, odd.
    """
    assert (N % 16) != 0, "Too many two factors"
    while N % 2 == 0:
        N = N // 2
    assert is_squarefree(N), "N not squarefree"
    return True


def group_data(N):
    """
    Computes necessary data associated to the group Gamma0(N).

    See the module description for groupdata for more.
    """
    check_level(N)
    congruence_subgroup = Gamma0(N)
    initial_vertices = [
        -1/2 + sqrt(3)/2 * i,
        1/2  + sqrt(3)/2 * i,
        infinity
    ]
    vertices = set()
    coset_reps = list(congruence_subgroup.coset_reps())
    # coset_reps = list(congruence_subgroup.farey_symbol().coset_reps())
    for gamma in coset_reps:
        for vertex in initial_vertices:
            vertices.add(gamma.acton(vertex))

    parabolic_vertices = set() # All parabolic vertices in fundamental domain, possibly
                               # with multiple corresponding to the same cusp.
    for vertex in vertices:
        if imag(vertex) == 0:
            parabolic_vertices.add(vertex)

    # reorder so that infinity is always the first cusp
    cusps = [infinity] + [QQ(cusp) for cusp in congruence_subgroup.cusps() if cusp != infinity]

    for cusp in cusps:
        if cusp == 0 or cusp == infinity:
            continue
        if cusp not in parabolic_vertices:
            raise ValueError("Cusp {} is not a pvertex in level {}".format(cusp, N))

    groupdata_dict = dict()
    groupdata_dict['coset_reps'] = coset_reps
    groupdata_dict['cusps'] = cusps
    groupdata_dict['pvertices'] = parabolic_vertices
    Udict = dict()
    Cdict = dict()
    for pvertex in parabolic_vertices:
        # (cusp, U_ell, sigma_ell)
        if pvertex == infinity:
            U_ell = matrix([[1, 0], [0, 1]])
            sigma = matrix([[1, 0], [0, 1]])
            Cdata = (infinity, infinity, U_ell, sigma)
            Cdict[infinity] = Cdata
            Udata = (infinity, infinity, U_ell, sigma)
            Udict[infinity] = Udata
        else:
            # This numerology is a pain. I briefly describe this on pg 45
            # of my notes, though it may be just as fast to directly compare
            # to Atkin-Lehner--Li.
            U_ell = congruence_subgroup.farey_symbol().reduce_to_cusp(pvertex)
            c = congruence_subgroup.reduce_cusp(pvertex)
            numer = c.numerator()
            denom = c.denominator()
            assert N % denom == 0
            width = N / denom
            g, s, t = xgcd(numer * width, -denom)
            assert g == 1
            sigma = matrix([[numer, t], [denom, s*width]])
            scaling = matrix([[sqrt(width), 0], [0, 1/sqrt(width)]])
            sigma = sigma * scaling
            if c not in Cdict:
                Cdata = (QQ(c), QQ(c), matrix([[1, 0], [0, 1]]), sigma)
                Cdict[QQ(c)] = Cdata
            Udata = (pvertex, QQ(c), U_ell, sigma)
            Udict[pvertex] = Udata
    groupdata_dict['Udata'] = Udict
    groupdata_dict['Cdata'] = Cdict
    groupdata_dict['level'] = N
    return groupdata_dict
