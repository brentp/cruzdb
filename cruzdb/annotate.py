from __future__ import print_function
import sys
import os
from cruzdb.models import Feature, ABase
from toolshed import reader, nopen

if sys.version_info[0] >= 3:
    basestring = str

def _annotate(args):
    print(args)
    try:
        return annotate(*args)
    except:
        print(args, file=sys.stderr)
        raise

def _split_chroms(fname):
    import tempfile
    t = tempfile.mktemp(dir="/tmp", suffix=".cruzdb")
    chroms = {}
    for d in reader(fname, header="ordered"):
        if not d['chrom'] in chroms:
            chroms[d['chrom']] = open(t + "." + d['chrom'], "w")
            print("\t".join(d.keys()), file=chroms[d['chrom']])
        print("\t".join(d.values()), file=chroms[d['chrom']])
    for k in chroms:
        chroms[k].close()
        chroms[k] = (chroms[k], chroms[k].name + ".anno")
    return chroms.items()

def annotate(g, fname, tables, feature_strand=False, in_memory=False,
        header=None, out=sys.stdout, _chrom=None, parallel=False):
    """
    annotate bed file in fname with tables.
    distances are integers for distance. and intron/exon/utr5 etc for gene-pred
    tables. if the annotation features have a strand, the distance reported is
    negative if the annotation feature is upstream of the feature in question
    if feature_strand is True, then the distance is negative if t
    """
    close = False
    if isinstance(out, basestring):
        out = nopen(out, "w")
        close = True


    if parallel:
        import multiprocessing
        import signal
        p = multiprocessing.Pool(initializer=lambda:
                                signal.signal(signal.SIGINT, signal.SIG_IGN))
        chroms = _split_chroms(fname)

        def write_result(fanno, written=[False]):
            for i, d in enumerate(reader(fanno, header="ordered")):
                if i == 0 and written[0] == False:
                    print("\t".join(d.keys()), file=out)
                    written[0] = True
                print("\t".join(x if x else "NA" for x in d.values()), file=out)
            os.unlink(fanno)
            os.unlink(fanno.replace(".anno", ""))

        for fchrom, (fout, fanno) in chroms:
            p.apply_async(annotate, args=(g.db, fout.name, tables, feature_strand, False, # True causes it to silently fail, resulting in nothing being writting in the output file
                                 header, fanno, fchrom),
                                 callback=write_result)
        p.close()
        p.join()
        return out.name

    if isinstance(g, basestring):
        from . import Genome
        g = Genome(g)
    if in_memory:
        from . intersecter import Intersecter
        intersecters = [] # 1 per table.
        for t in tables:
            q = getattr(g, t) if isinstance(t, basestring) else t
            if _chrom is not None:
                q = q.filter_by(chrom=_chrom)
            table_iter = q #page_query(q, g.session)
            intersecters.append(Intersecter(table_iter))

    elif isinstance(fname, basestring) and os.path.exists(fname) \
            and sum(1 for _ in nopen(fname)) > 25000:
        print("annotating many intervals, may be faster using in_memory=True", file=sys.stderr)
    if header is None:
        header = []
    extra_header = []
    for j, toks in enumerate(reader(fname, header=False)):
        if j == 0 and not header:
            if not (toks[1] + toks[2]).isdigit():
                header = toks
        if j == 0:
            for t in tables:
                annos = (getattr(g, t) if isinstance(t, basestring) else t).first().anno_cols
                h = t if isinstance(t, basestring) else t._table.name if hasattr(t, "_table") else t.first()._table.name
                extra_header += ["%s_%s" % (h, a) for a in annos]

            if 0 != len(header):
                if not header[0].startswith("#"):
                    header[0] = "#" + header[0]
                print("\t".join(header + extra_header), file=out)
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
            if len(objs) == 0:
                print("\t".join(toks + ["", "", ""]), file=out)
                continue

            gp = hasattr(objs[0], "exonStarts")
            names = [o.gene_name for o in objs]
            if feature_strand:
                strands = [-1 if f.is_upstream_of(o) else 1 for o in objs]
            else:
                strands = [-1 if o.is_upstream_of(f) else 1 for o in objs]

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

            # just take the first gene name if they are all the same
            if len(set(nd[0] for nd in name_dists)) == 1:
                toks.append(name_dists[0][0])
            else:
                toks.append(";".join(nd[0] for nd in name_dists))

            # iterate over the feat type, dist cols
            for i in range(1, len(name_dists[0])):

                toks.append(";".join(nd[i] for nd in name_dists))
        print("\t".join(toks), file=out)

    if close:
        out.close()
    return out.name
