from init import initialize_sql
from sqlalchemy import create_engine
import sys

from tests import test

def _open(filelike, mode='r'):
    if hasattr(filelike, 'read'): return filelike
    return open(filelike, mode)

class Genome(object):
    url = "mysql://%(user)s%(password)s@%(host)s/%(db)s"
    __tables = {}

    def __init__(self, db="", user="genome", host="genome-mysql.cse.ucsc.edu", password=""):
        if db.startswith(("sqlite://", "mysql://", "postgresql://")):
            self.url = db
            self.dburl = db
            self.user = self.host = self.password = ""
        else:
            self.db = db
            self.user = user
            self.host = host
            self.password = (":" + password) if password else ""
            self.dburl = self.url % dict(db=self.db, user=self.user,
                host=self.host, password=self.password)

        self.engine = create_engine(self.dburl)

        self.session, self.Base = initialize_sql(self.engine)
        self.models = __import__("models", globals(), locals(), ["Feature"], -1)

    def __del__(self):
        self.session.close_all()

    def mirror(self, tables, dest_url):
        from mirror import mirror
        mirror(self, tables, dest_url)

    def _map(self, table):
        # if the table hasn't been mapped, do so here.
        #if not table in self.__tables:
        if not table in self.Base.metadata.tables:
            # make a new class
            try:
                self.__tables[table] = type(table, (self.Base, getattr(self.models,
                    table)), {})
            except:
                self.__tables[table] = type(table, (self.Base,
                    self.models.Feature), {})
        return self.__tables[table]

    def __getattr__(self, table):
        self._map(table)
        mapped = self.session.query(self.__tables[table])
        mapped.table = lambda : self.table(table)
        mapped.orm = lambda : self._map(table)
        return mapped

    def table(self, table):
        self._map(table)
        return self.Base.metadata.tables[table]

    def bin_query(self, table, chrom, start, end):
        if isinstance(table, basestring):
            table = getattr(self, table)
        tbl = table.table()

        q = table
        q = q.filter(tbl.c.chrom == chrom)
        if hasattr(tbl.c, "bin"):
            bins = Genome.bins(start, end)
            q = q.filter(tbl.c.bin.in_(bins))
        if hasattr(tbl.c, "txStart"):
            return q.filter(tbl.c.txStart <= end, tbl.c.txEnd >= start)
        return q.filter(tbl.c.chromStart <= end, tbl.c.chromEnd >= start)

    def knearest(self, table, chrom_or_feat, start=None, end=None, k=1):

        if start is None:
            assert end is None
            chrom, start, end = chrom_or_feat.chrom, chrom_or_feat.start, chrom_or_feat.end
        else:
            chrom = chrom_or_feat

        qstart, qend = start, end
        res = self.bin_query(table, chrom, qstart, qend)
        change = 300
        while res.count() < k:
            qstart = max(0, qstart - change)
            qend += change
            change *= 2
            res = self.bin_query(table, chrom, qstart, qend)

        def dist(f):

            if start <= f.start <= end or start <= f.end <= end:
                d = 0
            else:
                d = min(abs(f.start - start),
                        abs(f.start - end),
                        abs(f.end - end),
                        abs(f.end - start))
            # add dist as an attribute to the feature
            f.dist = d
            return d


        # sort by dist and ...
        res = sorted(res, key=dist)
        if len(res) == k:
            return res
        ndist = res[k - 1].dist
        # include all features that are the same distance as the nth closest
        # feature (accounts for ties).
        while res[k - 1].dist == ndist and k < len(res):
            k = k + 1
        return res[:k]

    def sql(self, query):
        return self.engine.execute(query)

    @staticmethod
    def bins(start, end):
        if end - start < 536870912:
            offsets = [585, 73, 9, 1, 0]
        else:
            raise Exception("not implemented")
            offsets = [4681, 585, 73, 9, 1, 0]
        binFirstShift = 17
        binNextShift = 3

        start = start >> binFirstShift
        end = (end - 1)  >> binFirstShift

        bins = []
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
    g = Genome(db="hg18")#, host="localhost", user="")

    print g.cpgIslandExt[12].bed()
    print g.cpgIslandExt[12].bed('length', 'perCpg')

    import sys
    #sys.exit()

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
    query = g.refGene.filter(and_(g.table('refGene').c.txStart > 10000, g.table('refGene').c.txEnd < 40000))
    t = time.time()
    query.all()
    print time.time() - t

    query = g.refGene.filter(and_(g.table('refGene').c.txStart > 10000, g.table('refGene').c.txEnd < 40000))
    query = query.filter(g.table('refGene').c.bin.in_(Genome.bins(10000,
        40000)))

    t = time.time()
    query = g.bin_query(g.refGene, "chr1", 10000, 40000)

    query.all()
    print time.time() - t
    1/0


    g = Genome('hg19')
    t = time.time()
    q = g.snp135Common
    q = q.filter(q.table.c.bin.in_(Genome.bins(1000, 2000)))
    print q
    q.first()
    print time.time() - t

    #Genome.save_bed(query)
    1/0

    """
    for transcript in g.refGene:
        print transcript, transcript.sequence()[:100] + "..."
        if transcript.txEnd > 8000: break
    """

    kg = g.table('refGene')
    q = kg.select(kg.c.txStart < 5000)

    print list(g.session.execute(q))
