# from ordered_set import OrderedSet as oset
from scipy.special import binom
from sympy.combinatorics import Permutation

class Simplex:
    """

    """
    def __init__(faces):
        """
        faces: a dictionary with
        """
        self.faces = faces
        self.n = max(faces.keys())
        n = self.n
        for k in range(n):
            try:
                assert(binom(n,k) == len(faces[i]))
            except AssertionError as err:
                raise(err, "faces are not of the correct size!")
        return


    def intersect(other):
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
                new_sdict[i] |= self.sdict[i]

            return Simplex(new_sdict, o_n)

        # else, o_n >= s_n
        else:
            new_sdict = self.sdict.copy()

            # Union in all the parts from self
            for i in range(s_n+1):
                new_sdict[i] |= other.sdict[i]

            return Simplex(new_sdict, s_n)


    def __and__(other):
        """
        perform an operator override
        """
        return self.intersect(other)


    def union(other):
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

            return Simplex(new_sdict, s_n)

        # else, o_n >= s_n
        else:
            # If self is higher dimensional than other, then we'll
            # need a simplex with dimensions matching self's
            new_sdict = other.sdict.copy()

            # Union in all the parts from self
            for i in range(s_n+1):
                new_sdict[i] |= self.sdict[i]

            return Simplex(new_sdict, s_n)


    def __or__(other):
        return self.union(other)


class SimplicialComplex:
    """
    attributes: ordered sets of simplices?
    """
    def __init__(sdict, n):
        """
        initialize the simplex dictionary
        """
        self.sdict = sdict
        self.n = n
        return


    def intersect(other):
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
                new_sdict[i] |= self.sdict[i]

            return Simplex(new_sdict, o_n)

        # else, o_n >= s_n
        else:
            new_sdict = self.sdict.copy()

            # Union in all the parts from self
            for i in range(s_n+1):
                new_sdict[i] |= other.sdict[i]

            return Simplex(new_sdict, s_n)


    def __and__(other):
        """
        perform an operator override
        """
        return self.intersect(other)


    def union(other):
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

            return Simplex(new_sdict, s_n)

        # else, o_n >= s_n
        else:
            # If self is higher dimensional than other, then we'll
            # need a simplex with dimensions matching self's
            new_sdict = other.sdict.copy()

            # Union in all the parts from self
            for i in range(s_n+1):
                new_sdict[i] |= self.sdict[i]

            return Simplex(new_sdict, s_n)


    def __or__(other):
        return self.union(other)

def test_simplex():
    return
