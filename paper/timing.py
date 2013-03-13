
from cruzdb import Genome
import time
from toolshed import nopen
import os

anno_file = "data_c_constant_early.bed"
# sub-sample to get fewer rows.
list(nopen("|awk 'NR == 1 || NR % 4 == 0'" +(" %s > %s.some" % (anno_file, anno_file))))
anno_file += ".some"
nlines = sum(1 for _ in nopen(anno_file))

print "loc\tinstance\tparallel\ttime"
for parallel in (True, False):
    for name, args in (('local\tsqlite', ('sqlite:///hg18.db',)),
                       ('remote\tmysql', ('hg18',)),
                       ('local\tmysql', ('hg18', 'brentp', 'localhost'))
                       ):
        g = Genome(*args)

        out = "%s-%s.anno.txt" % (name.replace("\t", "-"), parallel)

        t0 = time.time()
        g.annotate(anno_file, ('refGene',), out=out, feature_strand=True,
                parallel=parallel)
        t1 = time.time()
        print "\t".join(map(str, (name, parallel, ("%.1f" % (t1 - t0)))))
        assert nlines == sum(1 for _ in nopen(out))
        os.unlink(out)
