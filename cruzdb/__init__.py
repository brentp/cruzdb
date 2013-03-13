import soup
import sys

class BigException(Exception): pass

from tests import test

def _open(filelike, mode='r'):
    if hasattr(filelike, 'read'): return filelike
    return open(filelike, mode)

class Genome(soup.Genome):
    url = "mysql://%(user)s%(password)s@%(host)s/%(db)s"

    def __init__(self, db="", user="genome", host="genome-mysql.cse.ucsc.edu",
            password="", engine=None):

        self.create_url(db, user, host, password)
        soup.Genome.__init__(self, self.dburl)
        self.session.autoflush = False

    def create_url(self, db="", user="genome", host="genome-mysql.cse.ucsc.edu",
        password=""):

        if db.startswith(("sqlite://", "mysql://", "postgresql://")):

            self.db = self.url = db
            self.dburl = db
            self.user = self.host = self.password = ""
        else:
            self.db = db
            if user == "genome" and host != "genome-mysql.cse.ucsc.edu":
                import getpass
                user = getpass.getuser()
            self.host = host
            self.user = user
            self.password = (":" + password) if password else ""
            self.dburl = self.url % dict(db=self.db, user=self.user,
                host=self.host, password=self.password)


    def mirror(self, tables, dest_url):
        from mirror import mirror
        return mirror(self, tables, dest_url)

    def dataframe(self, table, limit=None, offset=None):
        from pandas import DataFrame
        if isinstance(table, basestring):
            table = getattr(self, table)
        records = table._table.select()
        if not limit is None:
            records = records.limit(limit)
        if not offset is None:
            records = records.offset(offset)
        records = list(records.execute())
        cols = [c.name for c in table._table.columns]
        return DataFrame.from_records(records, columns=cols)

    @property
    def tables(self):
        return self._cache.keys()

    def load_file(self, fname, table=None, sep="\t", bins=False, indexes=None):
        """
        use some of the machinery in pandas to load a file into a table
        """
        convs = {"#chr": "chrom", "start": "txStart", "end": "txEnd", "chr":
                "chrom", "pos": "start", "POS": "start", "chromStart": "txStart",
                "chromEnd": "txEnd"}
        if table is None:
            import os.path as op
            table = op.basename(op.splitext(fname)[0]).replace(".", "_")
            print >>sys.stderr, "writing to:", table

        from pandas.io import sql
        import pandas as pa
        from toolshed import nopen

        needs_name = False
        for i, chunk in enumerate(pa.read_csv(nopen(fname), iterator=True,
            chunksize=100000, sep=sep, encoding="latin-1")):
            chunk.columns = [convs.get(k, k) for k in chunk.columns]
            if not "name" in chunk.columns:
                needs_name = True
                chunk['name'] = chunk.get('chrom', chunk[chunk.columns[0]])
            if bins:
                chunk['bin'] = 1
            if i == 0 and not table in self.tables:
                schema = sql.get_sqlite_schema(chunk, table)
                print schema
                self.engine.execute(schema)
            elif i == 0:
                print >>sys.stderr,\
                        """adding to existing table, you may want to drop first"""

            tbl = getattr(self, table)._table
            cols = chunk.columns
            data = list(dict(zip(cols, x)) for x in chunk.values)
            if needs_name:
                for d in data:
                    d['name'] = "%s:%s" % (d.get("chrom"), d.get("txStart", d.get("chromStart")))
            if bins:
                for d in data:
                    d['bin'] = max(Genome.bins(int(d["txStart"]), int(d["txEnd"])))
            self.engine.execute(tbl.insert(), data)
            self.session.commit()
            if i > 0:
                print >>sys.stderr, "writing row:", i * 100000
        if "txStart" in chunk.columns:
            if "chrom" in chunk.columns:
                ssql = """CREATE INDEX "%s.chrom_txStart" ON "%s" (chrom, txStart)""" % (table, table)
            else:
                ssql = """CREATE INDEX "%s.txStart" ON "%s" (txStart)""" % (table, table)

            self.engine.execute(ssql)
        for index in (indexes or []):
            ssql = """CREATE INDEX "%s.%s" ON "%s" (%s)""" % (table,
                                index, table, index)
            self.engine.execute(ssql)

        if bins:
            ssql = """CREATE INDEX "%s.chrom_bin" ON "%s" (chrom, bin)""" % (table, table)
            self.engine.execute(ssql)

        self.session.commit()

    @staticmethod
    def david_go(refseq_list, annot=('SP_PIR_KEYWORDS', 'GOTERM_BP_FAT',
                                        'GOTERM_CC_FAT', 'GOTERM_MF_FAT')):
        URL = "http://david.abcc.ncifcrf.gov/api.jsp?type=REFSEQ_MRNA&ids=%s&tool=term2term&annot="
        import webbrowser
        webbrowser.open(URL % ",".join(set(refseq_list)) + ",".join(annot))

    def bin_query(self, table, chrom, start, end):
        if isinstance(table, basestring):
            table = getattr(self, table)

        try:
            tbl = table._table
        except AttributeError:
            tbl = table.column_descriptions[0]['type']._table

        q = table.filter(tbl.c.chrom == chrom)

        if hasattr(tbl.c, "bin"):
            bins = Genome.bins(start, end)
            if len(bins) < 100:
                q = q.filter(tbl.c.bin.in_(bins))

        if hasattr(tbl.c, "txStart"):
            return q.filter(tbl.c.txStart <= end).filter(tbl.c.txEnd >= start)
        return q.filter(tbl.c.chromStart <= end).filter(tbl.c.chromEnd >= start)

    def upstream(self, table, chrom_or_feat, start=None, end=None, k=1):
        res = self.knearest(table, chrom_or_feat, start, end, k, "up")
        end = getattr(chrom_or_feat, "end", end)
        start = getattr(chrom_or_feat, "start", start)
        rev = getattr(chrom_or_feat, "strand", "+") == "-"
        if rev:
            return [x for x in res if x.end > start]
        else:
            return [x for x in res if x.start < end]

    def downstream(self, table, chrom_or_feat, start=None, end=None, k=1):
        res = self.knearest(table, chrom_or_feat, start, end, k, "down")
        end = getattr(chrom_or_feat, "end", end)
        start = getattr(chrom_or_feat, "start", start)
        rev = getattr(chrom_or_feat, "strand", "+") == "-"
        if rev:
            return [x for x in res if x.start < end]
        else:
            return [x for x in res if x.end > start]

    def knearest(self, table, chrom_or_feat, start=None, end=None, k=1,
            _direction=None):
        assert _direction in (None, "up", "down")

        # they sent in a feature
        if start is None:
            assert end is None
            chrom, start, end = chrom_or_feat.chrom, chrom_or_feat.start, chrom_or_feat.end

            # if the query is directional and the feature as a strand,
            # adjust...
            if _direction in ("up", "down") and getattr(chrom_or_feat,
                    "strand", None) == "-":
                _direction = "up" if _direction == "down" else "up"
        else:
            chrom = chrom_or_feat

        qstart, qend = long(start), long(end)
        res = self.bin_query(table, chrom, qstart, qend)

        i, change = 1, 350
        try:
            while res.count() < k:
                if _direction in (None, "up"):
                    if qstart == 0 and _direction == "up": break
                    qstart = max(0, qstart - change)
                if _direction in (None, "down"):
                    qend += change
                i += 1
                change *= (i + 5)
                res = self.bin_query(table, chrom, qstart, qend)
        except BigException:
            return []

        def dist(f):
            d = 0
            if start > f.end:
                d = start - f.end
            elif f.start > end:
                d = f.start - end
            # add dist as an attribute to the feature
            return d

        dists = sorted([(dist(f), f) for f in res])
        if len(dists) == 0:
            return []

        dists, res = zip(*dists)

        if len(res) == k:
            return res

        if k > len(res): # had to break because of end of chrom
            if k == 0: return []
            k = len(res)

        ndist = dists[k - 1]
        # include all features that are the same distance as the nth closest
        # feature (accounts for ties).
        while k < len(res) and dists[k] == ndist:
            k = k + 1
        return res[:k]

    def sql(self, query):
        return self.engine.execute(query)

    def annotate(self, fname, tables, feature_strand=False, in_memory=False,
            header=None, out=sys.stdout, parallel=False):
        from .annotate import annotate
        return annotate(self, fname, tables, feature_strand, in_memory, header=header,
                out=out, parallel=parallel)

    @staticmethod
    def bins(start, end):
        if end - start < 536870912:
            offsets = [585, 73, 9, 1]
        else:
            raise BigException
            offsets = [4681, 585, 73, 9, 1]
        binFirstShift = 17
        binNextShift = 3

        start = start >> binFirstShift
        end = (end - 1)  >> binFirstShift

        bins = [1]
        for offset in offsets:
            bins.extend(range(offset + start, offset + end + 1))
            start >>= binNextShift
            end >>= binNextShift
        return frozenset(bins)

    def __repr__(self):
        return "%s('%s')" % (self.__class__.__name__,
            self.url % dict(db=self.db, user=self.user, host=self.host,
            password="?" if self.password else ""))

    @classmethod
    def save_bed(cls, query, filename=sys.stdout):
        # write a bed12 file of the query.
        out = _open(filename, 'w')
        for o in query:
            out.write(o.bed() + '\n')


if __name__ == "__main__":
    #g = Genome(db="hg18", host="localhost", user="")
    import sys
    #g = Genome(db="hg18", host="localhost")
    g = Genome("sqlite:///hg18.db")

    #g.load_file('GSM882245.hg18.bed', bins=True)
    #g.load_file('http://rafalab.jhsph.edu/CGI/model-based-cpg-islands-hg18.txt',
    #        table='cpgRafaLab', bins=True)


    #1/0
    print g.cpgIslandExt[12].bed()
    print g.cpgIslandExt[12].bed('length', 'perCpg')

    #sys.exit()

    print "refGene"
    f = g.refGene[19]
    print f.bed12()
    f = g.refGene[19]
    print repr(f), f.cdsStart, f.cdsEnd
    print "exons", f.exons
    print "coding exons", f.coding_exons
    print "cds", f.cds

    print "introns", f.introns
    print "5'utr", f.utr5
    print "3'utr", f.utr3

    print f.browser_link
    #f.txEnd = f.txStart + 30
    #print list(f.blat())
    #print f.cds_sequence
    import time
    from sqlalchemy import and_
    query = g.refGene.filter(and_(g.refGene.txStart > 10000, g.refGene.txEnd < 40000))
    t = time.time()
    print query
    query.all()
    print time.time() - t

    query = g.refGene.filter(and_(g.refGene.txStart > 10000, g.refGene.txEnd < 40000))
    query = query.filter(g.refGene.bin.in_(Genome.bins(10000, 40000)))

    t = time.time()
    query = g.bin_query(g.refGene, "chr1", 10000, 40000)

    query.all()
    print time.time() - t


    g = Genome('hg19')
    t = time.time()
    q = g.snp135Common
    q = q.filter(q.bin.in_(Genome.bins(1000, 2000)))
    print q
    q.first()
    print time.time() - t

    Genome.save_bed(query)

    for transcript in g.refGene:
        print transcript, transcript.sequence()[:100] + "..."
        if transcript.txEnd > 8000: break

    kg = g.refGene._table
    q = kg.select(kg.c.txStart < 5000)

    print list(g.session.execute(q))
