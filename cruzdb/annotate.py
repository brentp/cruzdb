from cruzdb.models import Feature
import itertools
from toolshed import reader

def annotate(g, fname, tables):

    for i, toks in enumerate(reader(fname, header=False)):
        if i == 0 and not (toks[1] + toks[2]).isdigit():
            header = toks + ["%s_%s" % nv for nv in itertools.product(tables,
                                                        ("name", "distance"))]
            print "\t".join(header)
            continue
        f = Feature()
        f.chrom = toks[0]
        f.txStart = int(toks[1])
        f.txEnd = int(toks[2])

        for tbl in tables:
            objs = g.knearest(tbl, toks[0], int(toks[1]), int(toks[2]), k = 1)
            gp = objs[0].is_gene_pred
            toks.append(";".join(o.name for o in objs))
            toks.append(";".join(str(o.distance(f, features=gp)) for o in objs))
        print "\t".join(toks)

        if i > 200: break
