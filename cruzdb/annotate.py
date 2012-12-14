from cruzdb.models import Feature, ABase
import itertools
from toolshed import reader, nopen
import sys

def annotate(g, fname, tables, feature_strand=False, in_memory="auto",
        header=None, out=sys.stdout):
    """
    annotate bed file in fname with tables.
    distances are integers for distance. and intron/exon/utr5 etc for gene-pred
    tables. if the annotation features have a strand, the distance reported is
    negative if the annotation feature is upstream of the feature in question
    if feature_strand is True, then the distance is negative if t
    """

    if in_memory:
        from . intersecter import Intersecter
        from . mirror import page_query
        intersecters = []
        for t in tables:
            table_iter = page_query(getattr(g, t))
            intersecters.append(Intersecter(table_iter))
        #intersecters = [Intersecter(getattr(g, t)) for t in tables]

    elif isinstance(fname, basestring) and sum(1 for _ in nopen(fname)) > 2000:
        print >>sys.stderr, "annotating many intervals, may be faster using in_memory=True"
    if header is None:
        header = []
    extra_header = []

    for j, toks in enumerate(reader(fname, header=False)):
        if j == 0 and not header:
            if not (toks[1] + toks[2]).isdigit():
                header = toks

        if j == 0:
            for t in tables:
                annos = getattr(g, t).first().anno_cols
                extra_header += ["%s_%s" % (t, a) for a in annos]

            if not header[0].startswith("#"): header[0] = "#" + header[0]
            print >>out, "\t".join(header + extra_header)
            if header == toks: continue

        if not isinstance(toks, ABase):
            f = Feature()
            f.chrom = toks[0]
            f.txStart = int(toks[1])
            f.txEnd = int(toks[2])
            try:
                f.strand = toks[header.index('strand')]
            except ValueError:
                pass
        else:
            f = toks
            # for now, use the objects str to get the columns
            # might want to use getattr on the original cols

            toks = f.bed(*header).split("\t")

        sep = "^*^"
        for ti, tbl in enumerate(tables):
            if in_memory:
                objs = intersecters[ti].knearest(int(toks[1]), int(toks[2]), chrom=toks[0], k = 1)
            else:
                objs = g.knearest(tbl, toks[0], int(toks[1]), int(toks[2]), k=1)
            gp = hasattr(objs[0], "exonStarts")
            names = [o.gene_name for o in objs]
            if feature_strand:
                strands = [-1 if o.is_upstream_of(f) else 1 for o in objs]
            else:
                strands = [-1 if f.is_upstream_of(o) else 1 for o in objs]

            # dists can be a list of tuples where the 2nd item is something
            # like 'island' or 'shore'
            dists = [o.distance(f, features=gp) for o in objs]
            pure_dists = [d[0] if isinstance(d, (tuple, list)) else d for d in dists]

            # convert to negative if the feature is upstream of the query
            for i, s in enumerate(strands):
                if s == 1: continue
                if isinstance(pure_dists[i], basestring): continue
                pure_dists[i] *= -1

            for i, (pd, d) in enumerate(zip(pure_dists, dists)):
                if isinstance(d, tuple):
                    if len(d) > 1:
                        dists[i] = "%s%s%s" % (pd, sep, sep.join(d[1:]))
                    else:
                        dists[i] = pd
            # keep uniqe name, dist combinations (occurs because of
            # transcripts)
            name_dists = set(["%s%s%s" % (n, sep, d) \
                            for (n, d) in zip(names, dists)])
            name_dists = [nd.split(sep) for nd in name_dists]

            for i in range(len(name_dists[0])): # iterate over the dist, feature, name cols

                toks.append(";".join(nd[i] for nd in name_dists))
        print >>out, "\t".join(toks)
