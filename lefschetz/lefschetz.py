# For powerset generation
from itertools import combinations, chain

from scipy.special import binom
from sympy.combinatorics import Permutation


class SimplicialComplex:
    """
    attributes: ordered sets of simplices?
    """
    def __init__(self, sdict):
        """
        sdict:
            name: simplex dictionary
            type: dictionary with keys ints, values sets of ints representing vertices.
                {int : {verts}}

        """
        self.sdict = sdict
        self.n = max(sdict.keys())
        self.is_valid()
        return

    def __repr__(self):
        out_str = f"Simplicial Complex of max dimension {self.n}." + "\n---------------------------------------------\n" + "Subsimplices:\n"
        for i in range(self.n+1):
            out_str += f"    n={i}:"
            sorted_simps = sorted([sorted(list(simp)) for simp in list(self.sdict[i])])
            for i, simp in enumerate(sorted_simps):
                if i == 0:
                    out_str += " {" + str(simp)[1:-1] + "}\n"
                else:
                    # Trim out the list stuff
                    out_str += "         {" + str(simp)[1:-1] + "}\n"
            out_str += "\n"
        return out_str


    def is_valid(self):
        """
        check if self is a valid simplicial complex
        """
        bad_simps = {}
        # Check that all of the simplices are actually valid by
        for i in range(0, self.n+1, -1):
            for simplex in sdict[i]:
                is_valid, bad_simps = check_faces(simplex, bad_simps)
                try:
                    assert(is_valid)
                except AssertionError as err:
                    print(f"Invalid simplicial complex; not all faces of\
                            constituent simplices are included!\n\n Failed on\
                            simplex {simplex}.")
                    raise(err)
        return

    def check_faces(self, simplex, bad_simps):
        """
        Inputs:
        -------
        simplex: k-simplex (represented as a set of ints)
        """
        # Get the dimension of the simplex
        k = len(simplex)

        # Need to check all the way down to 0 simplices. k-simplices
        # aren't necessary, since this will only be called by
        # iterating over all k-simplices.
        for r in range(0, k, -1):
            r_combo_iterator = chain.from_iterable(combinations(simplex, r))
            for r_combo in r_subset_iterator:
                # chain.from_iterable spits out a tuple, not a set
                r_subset = set(r_combo)

                # Is this really a valid r simplex?
                if r_subset in self.sdict[r]:
                    pass
                elif r_subset in bad_simps:
                    return False, bad_simps
                else:
                    return False, bad_simps | {simplex, r_subset}
        return True, bad_simps

    def intersect(self, other):
        """
        iterate through the keys of the simplex dictionary and
        intersect the result
        """
        # self n, other simplex n
        s_n, o_n = self.n, other.n

        # new simplex dict
        new_sdict = {}
        if s_n > o_n:
            # If self is higher dimensional than other, then we'll
            # need a simplex with dimensions matching other's (because
            # the higher dimensional ones will intersect to empty)
            new_sdict = other.sdict.copy()

            # Intersect all the parts from self
            for i in range(o_n+1):
                new_sdict[i] &= self.sdict[i]

                # Delete empty sets
                if not new_sdict[i]:
                    del new_sdict[i]

            return SimplicialComplex(new_sdict)

        # o_n >= s_n
        else:
            new_sdict = self.sdict.copy()

            # Union in all the parts from self
            for i in range(s_n+1):
                new_sdict[i] &= other.sdict[i]

                # Delete empty sets
                if not new_sdict[i]:
                    del new_sdict[i]

            return SimplicialComplex(new_sdict)


    def __and__(self, other):
        """
        perform an operator override
        """
        return self.intersect(other)


    def union(self, other):
        """
        iterate through the simplex dictionaries
        """
        s_n, o_n = self.n, other.n
        # max_n = max(s_n, o_n)

        new_sdict = {}
        if s_n > o_n:
            # If self is higher dimensional than other, then we'll
            # need a simplex with dimensions matching self's
            new_sdict = self.sdict.copy()

            # Union in all the parts from other
            for i in range(o_n+1):
                new_sdict[i] |= other.sdict[i]

            return SimplicialComplex(new_sdict)

        # o_n >= s_n
        else:
            # If self is higher dimensional than other, then we'll
            # need a simplex with dimensions matching self's
            new_sdict = other.sdict.copy()

            # Union in all the parts from self
            for i in range(s_n+1):
                new_sdict[i] |= self.sdict[i]

            return SimplicialComplex(new_sdict)


    def __or__(self, other):
        return self.union(other)


    def subtract(self, other):
        """
        todo
        """
        return

# We need these things to work without having to go through the
# overhead of instantiating a class.
# TODO: refactor to define methods in terms of these, perhaps?
def is_valid(sdict):
    """
    check if self is a valid simplicial complex
    """
    bad_simps = {}
       # Check that all of the simplices are actually valid by
    for i in range(0, len(sdict.keys()), -1):
        for simplex in sdict[i]:
            is_valid, bad_simps = check_faces(sdict, simplex, bad_simps)
            if not is_valid:
                return False
    return True

def check_faces(sdict, simplex, bad_simps):
    """
    Inputs:
    -------
    simplex: k-simplex (represented as a set of ints)
    """
    # Get the dimension of the simplex
    k = len(simplex)

    # Need to check all the way down to 0 simplices. k-simplices
    # aren't necessary, since this will only be called by
    # iterating over all k-simplices.
    for r in range(0, k, -1):
        r_combo_iterator = chain.from_iterable(combinations(simplex, r))
        for r_combo in r_subset_iterator:
            # chain.from_iterable spits out a tuple, not a set
            r_subset = set(r_combo)

            # Is this really a valid r simplex?
            if r_subset in sdict[r]:
                pass
            elif r_subset in bad_simps:
                return False, bad_simps
            else:
                return False, bad_simps | {simplex, r_subset}
    return True, bad_simps

def apply_hom(f, K):
    """
    f: a dictionary associating verts in K to verts in the image of f

    Apply the simplicial map f over the simplex dict K
    """
    # Max value of n we can iterate to
    max_iter = len(K)+1
    out_sdict = {frozenset() for i in range(max_iter)}
    # Iterate through dimensions
    for n in range(max_iter):
        for knsimp in K[n]:
            # Image simplex
            imsimp = {}
            # Iterate through this n simplex in k
            for vert in knsimp:
                imsimp.add(f[vert])
            # Image of this n simplex need not be an n simplex
            out_sdict[len(imsimp)] |= frozenset(imsimp)
        # Currently in the middle of trying to write the logic to
        # merge all the images of the simplices (as determined by
        # verts) together
    # This needs testing!
    return out_sdict


def Hom_smap(K,L):
    """
    Homset of simplicial maps.

    K, L are simplicial complexes.

    generate all simplicial maps between the two simplicial complexes
    provided.
    """
    kverts, lverts = k.sdict[0], l.sdict[0]

    # Store previous vertex configurations that we've already tried in
    # a memo so that maybe we won't have n! runtime?
    memo = {}


    # Initialize output set of simplicial mpas
    all_smaps = {}

    # Current simplicial map (represented as a dictionary)
    smap = {}

    def smap_constructor(smap, kvset, lvset):
        """
        kvset, lvset --- same as kverts, lverts. Just using a
        different name b/c I forgot how global vars work in python and
        want to avoid naming conflicts.
        """
        # Base case: all points in the domain have been assigned.
        if not(kverts):
            return smap
        for bad_map in memo:
            if dict_subset(bad_map, smap):
                return False

        else:
            for lv in lvset:
                new_smap = smap.copy()

                # pop a vertex out
                new_kvset = kvset.copy()
                this_k = new_kvset.pop()

                # Randomly assign an image
                new_smap[this_k] = lv

                new_complex = apply_hom(new_smap, kvset)
                if is_valid(new_complex):
                    # TODO: Maybe figure out a way to lazy-evaluate
                    # these things. Perhaps an iterator could be
                    # useful?

                    return smap_constructor(smap_copy, new_kvset, lvset)
                else:
                    memo.add(smap)
                    return False

    return

def dict_subset(d1, d2):
    """
    check if one of the input dictionaries is a subset of the other.

    in particular, if d1 <= d2
    """
    return all(d1.get(k, object()) == v for (k, v) in d2.items())




def standard_simp(verts):
    """
    For recursively constructing a standard n-simplex from a list of verts.
    Basically, so that I can pass in something like [2,3,4], and get all of the
    things on those vertices without having to go through the construction
    explicitly.

    returns a simplex dictionary
    """
    k = len(verts)

    # Initialize the dictionary
    sdict = {}

    for r in range(1,k+1):
        sub_faces = {frozenset(r_face) for r_face in combinations(verts, r)}
        # 0-simplices are points, not the empty set.
        sdict[r-1] = sub_faces
    return sdict


def test_simplex():
    sdict1 = standard_simp({1,2,3})
    sdict2 = standard_simp({2,3,4})
    sdict3 = standard_simp({5,6,7,8})

    scomp1 = SimplicialComplex(sdict1)
    scomp2 = SimplicialComplex(sdict2)
    scomp3 = SimplicialComplex(sdict3)
    return scomp1, scomp2, scomp3

scomp1, scomp2, scomp3 = test_simplex()
