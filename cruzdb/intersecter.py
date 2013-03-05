import operator
import collections

class Feature(object):
    """\
    Basic feature, with required integer start and end properties.
    Also accpets optional strand as +1 or -1 (used for up/downstream queries),
    a name

    >>> from intersecter import Feature

    >>> f1 = Feature(23, 36)
    >>> f2 = Feature(34, 48, strand=-1)
    >>> f2
    Feature(34, 48, strand=-1)

    """
    __slots__ = ("start", "end", "strand", "chrom")

    def __init__(self, start, end, strand=0, chrom=None):
        assert start <= end, "start must be less than end"
        self.start  = start
        self.end   = end
        self.strand = strand
        self.chrom  = chrom

    def __repr__(self):
        fstr = "Feature(%d, %d" % (self.start, self.end)
        if self.chrom is not None:
            fstr += ", chrom=%s" % self.chrom
        if self.strand != 0:
            fstr += ", strand=%d" % self.strand
        fstr += ")"
        return fstr

def binsearch_left_start(intervals, x, lo, hi):
    while lo < hi:
        mid = (lo + hi)//2
        f = intervals[mid]
        if f.start < x: lo = mid + 1
        else: hi = mid
    return lo

# like python's bisect_right find the _highest_ index where the value x 
# could be inserted to maintain order in the list intervals
def binsearch_right_end(intervals, x, lo, hi):
    while lo < hi:
        mid = (lo + hi)/2
        f = intervals[mid]
        if x < f.start: hi = mid
        else: lo = mid + 1
    return lo

class Intersecter(object):
    """\
    Data structure for performing intersect and neighbor queries on a
    set of intervals. Algorithm uses simple binary search along with
    knowledge of the longest interval to perform efficient queries.

    Usage
    =====
    >>> from intersecter import Intersecter, Feature

    Add intervals, the only requirement is that the interval have integer
    start and end attributes. Optional arguments are strand, and chrom.

    >>> f = Feature(1, 22, strand=-1)
    >>> f
    Feature(1, 22, strand=-1)

    >>> features = [
    ...            Feature(0, 10, -1),
    ...            Feature(3, 7, 1),
    ...            Feature(3, 40, -1),
    ...            Feature(13, 50, 1)
    ... ]

    >>> intersecter = Intersecter(features)

    Queries
    -------

    find
    ++++

    >>> intersecter.find(2, 5)
    [Feature(0, 10, strand=-1), Feature(3, 7, strand=1), Feature(3, 40, strand=-1)]
    >>> intersecter.find(11, 100)
    [Feature(3, 40, strand=-1), Feature(13, 50, strand=1)]
    >>> intersecter.find(100, 200)
    []

    left/right
    ++++++++++
    the left method finds features that are strictly to the left of
    the query feature. overlapping features are not considered:

    >>> intersecter.left(Feature(0, 1))
    []
    >>> intersecter.left(Feature(11, 12))
    [Feature(0, 10, strand=-1)]


    up/downstream
    +++++++++++++
    up/downstream method behave exactly like left/right, except that
    the direction is determined by the strand of the query feature. 
    If the strand is 1, then upstream is left, downstream is right.

    If the strand is -1, then upstream is right, downstream is left.
    >>> intersecter.upstream(Feature(11, 12, strand=1))
    [Feature(0, 10, strand=-1)]
    >>> intersecter.upstream(Feature(11, 12, strand=-1))
    [Feature(13, 50, strand=1)]

    all of these method take an argument 'n' for the number of results desired.
    >>> intersecter.upstream(Feature(1, 2, strand=-1), n=3)
    [Feature(3, 7, strand=1), Feature(3, 40, strand=-1), Feature(13, 50, strand=1)]

    nearest neighbors
    +++++++++++++++++
    >>> intersecter.knearest(Feature(1, 2))
    [Feature(0, 10, strand=-1)]

    >>> intersecter.knearest(Feature(1, 2), k=2)
    [Feature(0, 10, strand=-1), Feature(3, 40, strand=-1), Feature(3, 7, strand=1)]

    """

    # since intervals are sorted by start, also have to know the max_len (see find)
    # cdef int max_len
    # if an item is added, the list must be resorted.

    # ---- Basic API --------------------------------------------------

    def __init__(self, intervals):
        self.intervals = collections.defaultdict(list)
        self.max_len = {}

        for i, iv in enumerate(intervals):
            self.intervals[iv.chrom].append(iv)
        for chrom in self.intervals:
            self.intervals[chrom].sort(key=operator.attrgetter('start'))
            self.max_len[chrom] = max(1, max([i.end - i.start for i in self.intervals[chrom]]))


    def find(self, start, end, chrom=None):
        """Return a object of all stored intervals intersecting between (start, end) inclusive."""
        intervals = self.intervals[chrom]
        ilen = len(intervals)
        # NOTE: we only search for starts, since any feature that starts within max_len of
        # the query could overlap, we must subtract max_len from the start to get the needed
        # search space. everything else proceeds like a binary search.
        # (but add distance calc for candidates).
        if not chrom in self.max_len: return []
        ileft  = binsearch_left_start(intervals, start - self.max_len[chrom], 0, ilen)
        iright = binsearch_right_end(intervals, end, ileft, ilen)
        query = Feature(start, end)
        # we have to check the distance to make sure we didnt pick up anything 
        # that started within max_len, but wasnt as long as max_len
        return [f for f in intervals[ileft:iright] if distance(f, query) == 0]

    def left(self, f, n=1):
        """return the nearest n features strictly to the left of a Feature f.
        Overlapping features are not considered as to the left.

        f: a Feature object
        n: the number of features to return
        """
        intervals = self.intervals[f.chrom]
        if intervals == []: return []

        iright = binsearch_left_start(intervals, f.start, 0 , len(intervals)) + 1
        ileft  = binsearch_left_start(intervals, f.start - self.max_len[f.chrom] - 1, 0, 0)

        results = sorted((distance(other, f), other) for other in intervals[ileft:iright] if other.end < f.start and distance(f, other) != 0)
        if len(results) == n:
            return [r[1] for r in results]

        # have to do some extra work here since intervals are sorted
        # by starts, and we dont know which end may be around...
        # in this case, we got some extras, just return as many as
        # needed once we see a gap in distances.
        for i in range(n, len(results)):
            if results[i - 1][0] != results[i][0]:
                return [r[1] for r in results[:i]]

        if ileft == 0:
            return [r[1] for r in results]

        # here, didn't get enough, so move left and try again. 
        1/0

    def right(self, f, n=1):
        """return the nearest n features strictly to the right of a Feature f.
        Overlapping features are not considered as to the right.

        f: a Feature object
        n: the number of features to return
        """
        intervals = self.intervals[f.chrom]
        ilen = len(intervals)
        iright = binsearch_right_end(intervals, f.end, 0, ilen)
        results = []

        while iright < ilen:
            i = len(results)
            if i > n:
                if distance(f, results[i - 1]) != distance(f, results[i - 2]):
                    return results[:i - 1]
            other = intervals[iright]
            iright += 1
            if distance(other, f) == 0: continue
            results.append(other)
        return results


    def upstream(self, f, n=1):
        """find n upstream features where upstream is determined by
        the strand of the query Feature f
        Overlapping features are not considered.

        f: a Feature object
        n: the number of features to return
        """
        if f.strand == -1:
            return self.right(f, n)
        return self.left(f, n)


    def downstream(self, f, n=1):
        """find n downstream features where downstream is determined by
        the strand of the query Feature f
        Overlapping features are not considered.

        f: a Feature object
        n: the number of features to return
        """
        if f.strand == -1:
            return self.left(f, n)
        return self.right(f, n)

    def knearest(self, f_or_start, end=None, chrom=None, k=1):
        """return the n nearest neighbors to the given feature
        f: a Feature object
        k: the number of features to return
        """


        if end is not None:
            f = Feature(f_or_start, end, chrom=chrom)
        else:
            f = f_or_start

        DIST = 2000
        feats = filter_feats(self.find(f.start - DIST, f.end + DIST, chrom=f.chrom), f, k)
        if len(feats) >= k:
            return feats

        nfeats = k - len(feats)
        fleft = Feature(f.start - DIST, f.start, chrom=f.chrom)
        feats.extend(self.left(fleft, n=nfeats))

        fright = Feature(f.end, f.end + DIST, chrom=f.chrom)
        feats.extend(self.right(fright, n=nfeats))
        return filter_feats(feats, f, k)


def distance(f1, f2):
    """\
    Distance between 2 features. The integer result is always positive or zero.
    If the features overlap or touch, it is zero.
    >>> from intersecter import Feature, distance
    >>> distance(Feature(1, 2), Feature(12, 13))
    10
    >>> distance(Feature(1, 2), Feature(2, 3))
    0
    >>> distance(Feature(1, 100), Feature(20, 30))
    0

    """
    if f1.end < f2.start: return f2.start - f1.end
    if f2.end < f1.start: return f1.start - f2.end
    return 0

def filter_feats(intervals, f, k):
    feats = sorted((distance(f, iv), iv) for iv in intervals if iv is not None)
    kk = k
    while kk < len(feats) and feats[k - 1][0] == feats[kk][0]:
        kk += 1
    return [f[1] for f in feats[:kk]]

if __name__ == "__main__":
    import doctest
    doctest.testmod()
