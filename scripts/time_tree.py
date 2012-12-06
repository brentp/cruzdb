import sys
sys.path.extend([".", "scripts", "cruzdb"])
import operator
import random
from interval_tree import Interval, IntervalTree
from intersecter import Intersecter
import time


N=500000
TRIES = 10
START, STOP = 0, 20000

def rands(n=N, len_range=(200, 16000), start_max=STOP):

    def rand():
        start = random.randint(1, start_max)
        return Interval(start, start + random.randint(*len_range))

    return [rand() for i in xrange(n)]


def brute_force_find(intervals, start, end):
    return [i for i in intervals if i.end >= start and i.start <= end]

def search(tree, start, end, tries):

    t0 = time.time()
    lens = []
    if isinstance(tree, list):
        for i in range(tries):
            res = brute_force_find(tree, start, end)
            res.sort(key=operator.attrgetter('start'))
            lens.append("%i:%s" % (len(res), [x.start for x in res[-1:]]))
            #lens.append(len(res))
    else:
        for i in range(tries):
            res = tree.find(start, end)
            res.sort(key=operator.attrgetter('start'))
            lens.append("%i:%s" % (len(res), [x.start for x in res[-1:]]))
            #lens.append(len(res))
    t1 = time.time()
    return res, t1 - t0, lens


start_max = STOP * 3
while True:
    intervals = rands(N, start_max = start_max)
    t0 = time.time()
    tree = IntervalTree(intervals)
    t1 = time.time()
    print "time to build IntervalTree with %i intervals: %.3f" % (N, t1 - t0)
    t0 = time.time()
    ints = Intersecter(intervals)
    t1 = time.time()
    print "time to build Intersector with %i intervals: %.3f" % (N, t1 - t0)

    found, t, tree_lens = search(tree, START, STOP, TRIES)
    print "time to search tree %i times: %.3f. found %i intervals" % (TRIES, t, len(found))

    found, t, brute_lens = search(intervals, START, STOP, TRIES)
    print "time to search brute %i times: %.3f. found %i intervals" % (TRIES, t, len(found))

    found, t, inter_lens = search(ints, START, STOP, TRIES)
    print "time to search intersecter %i times: %.3f. found %i intervals" % (TRIES, t, len(found))

    for tl, bl, il in zip(tree_lens, brute_lens, inter_lens):
        assert tl == bl == il, (tl, bl, il)
    print
    #start_max *= 2
