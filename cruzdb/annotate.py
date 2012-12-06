from cruzdb.models import Feature
import itertools
from toolshed import reader

def annotate(g, fname, tables, feature_strand=False):
    """
    annotate bed file in fname with tables.
    distances are integers for distance. and intron/exon/utr5 etc for gene-pred
    tables. if the annotation features have a strand, the distance reported is
    negative if the annotation feature is upstream of the feature in question
    if feature_strand is True, then the distance is negative if t
    """

    for j, toks in enumerate(reader(fname, header=False)):
        if j == 0 and not (toks[1] + toks[2]).isdigit():

            header = toks + ["%s_%s" % nv for nv in itertools.product(tables,
                                                        ("name", "distance"))]
            print "\t".join(header)
            continue
        f = Feature()
        f.chrom = toks[0]
        f.txStart = int(toks[1])
        f.txEnd = int(toks[2])
        f.strand = toks[header.index('strand')]

        sep = "^*^"
        for tbl in tables:
            objs = g.knearest(tbl, toks[0], int(toks[1]), int(toks[2]), k = 1)
            gp = objs[0].is_gene_pred
            names = [o.gene_name for o in objs]
            if feature_strand:
                strands = [-1 if o.is_upstream_of(f) else 1 for o in objs]
            else:
                strands = [-1 if f.is_upstream_of(o) else 1 for o in objs]

            # dists can be a list of tuples where the 2nd item is something
            # like 'island' or 'shore'
            dists = [o.distance(f, features=gp) for o in objs]
            pure_dists = [d[0] if isinstance(d, tuple) else d for d in dists]

            # convert to negative if the feature is upstream of the query
            for i, s in enumerate(strands):
                if s == 1: continue
                if isinstance(pure_dists[i], basestring): continue
                pure_dists[i] *= -1

            for i, (pd, d) in enumerate(zip(pure_dists, dists)):
                if isinstance(d, tuple):
                    if d[1] not in ("", None):
                        dists[i] = "%s/%s" % (pd, d[1])
                    else:
                        dists[i] = pd
            # keep uniqe name, dist combinations (occurs because of
            # transcripts)
            name_dists = set(["%s%s%s" % (n, sep, d) \
                            for (n, d) in zip(names, dists)])
            name_dists = [nd.split(sep) for nd in name_dists]

            toks.append(";".join(nd[0] for nd in name_dists))
            toks.append(";".join(nd[1] for nd in name_dists))
        print "\t".join(toks)

        if j > 200: break
