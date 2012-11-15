"""
"testing" of bin queries and nearest queries
take a long time to run so not part of standard test suite
"""

import time
from cruzdb import Genome
from random import randrange, seed

#g = Genome('hg18', host='localhost', user='brentp')

#g.mirror(['refGene'], "sqlite:////tmp/u.db")


g = Genome('sqlite:////tmp/u.db')



# if we choose a huge distance all should have a distance of 0
#assert all(k.dist == 0  for k in  g.knearest("refGene", "chr1", 1234, 9915555, k=3))


print g.upstream("refGene", "chr1", 9444, 9555, k=6)

last = g.refGene.order_by(-g.refGene.table().c.txStart)[0]
print last
last.txStart = 1000 + last.txEnd
last.txEnd = last.txStart + 100
last.strand = "-"
print last
print g.upstream("refGene", last, k=6)

1/0


seed(1)

istart = 12345
iend = 386539

qall = list(g.refGene.all())

#while True:
for iend in (randrange(istart, 65555555) for i in range(100)):


    t = time.time()
    q = g.bin_query('refGene', 'chr1', istart, iend)
    a = list(q)
    print len(a)
    print time.time() - t

    #"""
    t = time.time()
    refGene = g.refGene


    rg = refGene.table()
    q = g.session.query(rg).filter(rg.c.chrom == "chr1", rg.c.txStart
            <= iend, rg.c.txEnd >= istart)
    q = refGene.filter(rg.c.chrom == "chr1", rg.c.txStart
            <= iend, rg.c.txEnd >= istart)
    b = list(q)
    print len(b)

    print time.time() - t
    #"""

    t = time.time()
    b = [r for r in g.refGene.all() if r.chrom == "chr1" and r.txStart <= iend and r.txEnd >= istart]
    #print(len(qall))
    #b = [r for r in qall if (r.chrom == "chr1" and r.txStart <= iend and r.txEnd >= istart)]
    print time.time() - t

    assert len(a) == len(b), (len(a), len(b), iend)
    print
