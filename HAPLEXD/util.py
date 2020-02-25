import pandas as pd
import numpy as np
import functools


def read_haplotypes(filename):
    haplotypeDF= pd.read_csv(filename, header="infer")
    haplotypesDF=haplotypesDF.transpose()
    return haplotypesDF.to_numpy()

def tractize(haplotype):
    """
    Given a haplotype sequence, returns the tractized version
    """
    tractized=[]
    for i in range(len(haplotype)):
        tractized.append((2*i)+haplotype[i])
    return tractized

@functools.total_ordering
class UniqueEndChar (object):
    """ 
    Termination character. 
    This class helps create a unique termination character
    for each haplotype given in the training sample.
     """
    def __init__ (self,hap_id):
        self.id = hap_id
    def __str__ (self):
        return '$'
    def __lt__ (self, other):
        return False


class Path (object):
    """ A path in a suffix tree. """
    e = 0
    inf = -1
    def __init__ (self, S, start = 0, end = None):
        if end is None:
            end = len(S)

        _end = end if end >= 0 else self.e
        assert 0 <= start <= _end <= len (S), "Path: 0 <= %d <= %d <= %d" % (start, end, len (S))

        self.S     = S
        self.start = start
        self._end  = end

    @property
    def end (self):
        """
        This needs its own calculation to be set as a property because the end can sometimes be open.
        """
        if self._end!= self.inf:
            return self._end
        else:
            return self.e
        
    @end.setter
    def end (self, value):
        self._end = value

    @classmethod
    def construct (cls, path_info):
        path_info = tuple(path_info)
        return Path (path_info, 0, len (path_info))

    def __str__ (self):
        return ' '.join (map (str, self.S[self.start:self.end]))

    def __len__ (self):
        return self.end - self.start

    def __lt__ (self, other):
        return self.S[self.start:self.end] < other.S[other.start:other.end]

    def __getitem__ (self, key):
        length = len(self)
        if isinstance (key, slice):
            start, stop, step = key.indices (length)
            if step != 1:
                raise TypeError
            return self.S[self.start + start:self.start + stop]
        if isinstance (key, int):
            if key < 0:
                key = len (self) + key
            if key >= length:
                raise IndexError ("key = %d, len () = %d" % (key, length))
            return self.S[self.start + key]
        raise TypeError ()
