# For powerset generation
from copy import deepcopy

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
        # assert(self.is_valid())
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
        bad_simps = set()
        # Check that all of the simplices are actually valid by
        for i in range(self.n, -1, -1):
            for simplex in self.sdict[i]:
                valid, bad_simps = self.check_faces(simplex, bad_simps, i)
                if not valid:
                    print("Invalid simplicial complex; not all faces of constituent simplices are included!\n\nFailed on simplex " + str(simplex)[10:-1] + " in the following complex:")
                    print(self)
                    print("Complex is missing ", str(bad_simps)[10:-1], "\n\n\n\n")
                    return False
        return True

    def check_faces(self, simplex, bad_simps, k):
        """
        Inputs:
        -------
        simplex: k-simplex (represented as a set of ints)
        """
        # Need to check all the way down to 0 simplices. k-simplices
        # aren't necessary, since this will only be called by
        # iterating over all k-simplices.

        # TODO --- bad_simps sort of isn't necessary here, is it?
        for r in range(k-1, -1, -1):
            r_combo_iterator = chain.from_iterable(combinations(simplex, r))
            for r_combo in r_combo_iterator:
                # chain.from_iterable spits out a tuple, not a set
                r_subset = frozenset({r_combo})
                i = len(r_subset)-1

                # Is this really a valid r simplex?
                if r_subset in self.sdict[i]:
                    pass
                elif r_subset in bad_simps:
                    return False, bad_simps
                else:
                    # Or the r_subset into bad_simps
                    return False, bad_simps | {r_subset}
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
            new_sdict = deepcopy(other.sdict)

            # Intersect all the parts from self
            for i in range(o_n+1):
                new_sdict[i] &= self.sdict[i]

                # Delete empty sets
                if not new_sdict[i]:
                    del new_sdict[i]

            return SimplicialComplex(new_sdict)

        # o_n >= s_n
        else:
            new_sdict = deepcopy(self.sdict)

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
            new_sdict = deepcopy(self.sdict)

            # Union in all the parts from other
            for i in range(o_n+1):
                new_sdict[i] |= other.sdict[i]

            return SimplicialComplex(new_sdict)

        # o_n >= s_n
        else:
            # If self is higher dimensional than other, then we'll
            # need a simplex with dimensions matching self's
            new_sdict = deepcopy(other.sdict)

            # Union in all the parts from self
            for i in range(s_n+1):
                new_sdict[i] |= self.sdict[i]

            return SimplicialComplex(new_sdict)


    def __or__(self, other):
        return self.union(other)


    def subtract(self, other):
        """
        TODO
        """
        return

# We need these things to work without having to go through the
# overhead of instantiating a class.
# TODO: refactor

def valid_sdict(sdict):
    """
    check if self is a valid simplicial complex
    """
    bad_simps = set()
    # Check that all of the simplices are actually valid by
    for i in range(len(sdict)-1, -1, -1):
        for simplex in sdict[i]:
            valid, bad_simps = sdict_check_faces(sdict, simplex, bad_simps, i)
            if not valid:
                # print("Invalid simplicial complex; not all faces of constituent simplices are included!\n\nFailed on simplex " + str(simplex)[10:-1] + " in the following complex:")
                # print(sdict)
                # print("Complex is missing ", str(bad_simps)[10:-1], "\n\n\n\n")
                return False
    return True

def sdict_check_faces(sdict, simplex, bad_simps, k):
    """
    Inputs:
    -------
    simplex: k-simplex (represented as a set of ints)
    """
    # Need to check all the way down to 0 simplices. k-simplices
    # aren't necessary, since this will only be called by
    # iterating over all k-simplices.

    # TODO --- bad_simps sort of isn't necessary here, is it?
    for r in range(k-1, -1, -1):
        r_combo_iterator = chain.from_iterable(combinations(simplex, r))
        for r_combo in r_combo_iterator:
            # chain.from_iterable spits out a tuple, not a set
            r_subset = frozenset({r_combo})
            i = len(r_subset)-1

            # Is this really a valid r simplex?
            if r_subset in sdict[i]:
                pass
            elif r_subset in bad_simps:
                return False, bad_simps
            else:
                # Or the r_subset into bad_simps
                return False, bad_simps | {r_subset}
    return True, bad_simps



def apply_hom(f, K):
    """
    f: a dictionary associating verts in K to verts in the image of f

    Apply the simplicial map f over the simplex dict K
    """
    # Max value of n we can iterate to
    max_iter = len(K)-1
    # print(max_iter)
    out_sdict = {i:set() for i in range(max_iter+1)}
    # Iterate through dimensions
    for n in range(max_iter-1, -1, -1):
        for knsimp in K[n]:
            # Image simplex
            imsimp = set()
            # Iterate through this n simplex in k
            for vert in knsimp:
                if vert in f:
                    imsimp.add(f[vert])

            # Image of this n simplex need not be an n simplex
            try:
                out_sdict[len(imsimp)] |= {frozenset(imsimp)}
            except KeyError as err:
                print("out_sdict: ", out_sdict)
                print("imsimp: ", imsimp)
                print("key attempted: ", len(imsimp)-1)
                raise(err)

    # Filter out empty entries
    out_sdict = {k:v for (k,v) in out_sdict.items() if v}
    return out_sdict

def smap_constructor(memo, all_smaps, smap, kverts, lverts, ksdict):
    # kverts, lverts = K.sdict[0], L.sdict[0]

    if not kverts:
        return memo, all_smaps
    for bad_map in memo:
        if frozendict_subset(bad_map, smap):
            return False

    # Iterate all k vertices --- this actually double-counts still. How to fix this?
    for kv in kverts:
        new_kverts = deepcopy(kverts)
        new_kverts.remove(kv)
        # Iterate
        for lv in lverts:
            new_smap = deepcopy(smap)
            new_smap[kv] = lv
            if valid_sdict(apply_hom(new_smap, ksdict)):
                next_level = smap_constructor(memo, all_smaps, new_smap, new_kverts, lverts, ksdict)
                if next_level:
                    print(next_level)
                    memo, all_smaps = next_level
                    # Fricken python not having a immutable dictionary type... smh
                    all_smaps.add(frozenset(new_smap.items()))
            else:
                memo.add(frozenset(new_smap.items()))
                return False

    return memo, all_smaps


def Hom(K,L):
    """
    Homset of simplicial maps.

    K, L are simplicial complexes.

    generate all simplicial maps between the two simplicial complexes
    provided.
    """
    kverts, ksdict = K.sdict[0], K.sdict
    lverts = L.sdict[0]
    # We don't need the memo
    _, all_smaps = smap_constructor(set(), set(), dict(), kverts, lverts, ksdict)
    return all_smaps

def dict_subset(d1, d2):
    """
    check if one of the input dictionaries is a subset of the other.

    in particular, if d1 <= d2
    """
    d1 = dict(d1)
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

    s1dict = {0 : {frozenset({1})},
              1 : {frozenset({1,2})},
              2 : {frozenset({1,2,3})}}

    s2dict = {0 : {frozenset({1}), frozenset({2})},
              1 : {frozenset({1,2}), frozenset({2,3}), frozenset({1,3})},
              2 : {frozenset({1,2,3})}}

    s3dict = {0 : {frozenset({1}), frozenset({2}), frozenset({3}), frozenset({4})},
              1 : {frozenset({1,2}), frozenset({2,3}), frozenset({1,3})}}


    s1 = SimplicialComplex(s1dict)
    s2 = SimplicialComplex(s2dict)
    s3 = SimplicialComplex(s3dict)


    assert(scomp1.is_valid())
    assert(not s1.is_valid())
    assert((scomp1 & scomp2).is_valid())
    assert(not s2.is_valid())
    assert(s3.is_valid())

    return scomp1, scomp2, scomp3, s1

scomp1, scomp2, scomp3, s1 = test_simplex()

f = {1:1, 2:1, 3:2}
g = {1:1, 2:3, 3:2}
kp = apply_hom(f, scomp1.sdict)
lp = apply_hom(g, scomp1.sdict)

aut_scomp1 = Hom(scomp1, scomp1)
