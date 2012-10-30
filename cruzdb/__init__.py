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
        self.db = db
        self.user = user
        self.host = host
        self.password = (":" + password) if password else ""
        self.dburl = self.url % dict(db=self.db, user=self.user,
            host=self.host, password=self.password)
        self.engine = create_engine(self.dburl)

        self.session, self.Base = initialize_sql(self.engine)
        self.models = __import__("models", globals(), locals(), ["Feature"], -1)

    def _map(self, table):
        # if the table hasn't been mapped, do so here.
        #if not table in self.models.Base.metadata.tables:
        if not table in self.__tables:
            # make a new class
            try:
                self.__tables[table] = type(table, (self.Base, getattr(self.models, table)), {})
            except:
                self.__tables[table] = type(table, (self.Base,
                    self.models.Feature), {})

    def __getattr__(self, table):
        self._map(table)
        return self.session.query(self.__tables[table])

    def table(self, table):
        self._map(table)
        return self.Base.metadata.tables[table]

    def sql(self, query):
        return self.engine.execute(query)

    @staticmethod
    def bins(start, end):
        bins = [1]
        bins.extend(xrange(1 + (start>>26), 1 + ((end-1)>>26)+1))
        bins.extend(xrange(9 + (start>>23), 9 + ((end-1)>>23)+1))
        bins.extend(xrange(73 + (start>>20), 73 + ((end-1)>>20)+1))
        bins.extend(xrange(585 + (start>>17), 585 + ((end-1)>>17)+1))
        return frozenset(bins)


    @staticmethod
    def bin_query(start, end, query):
        return query.filter.c.bin.in_(Genome.bins(start, end))

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
    #print f.cds_sequence
    import time
    from sqlalchemy import and_
    query = g.knownGene.filter(and_(g.table('knownGene').c.txStart > 10000, g.table('knownGene').c.txEnd < 40000))
    t = time.time()
    query.all()
    print time.time() - t

    query = query.filter(g.table('knownGene').c.bin.in_(Genome.bins(10000,
        40000)))

    t = time.time()
    query = Genome.bin_query(10000, 40000, g.table('knownGene'))
    query.all()
    print time.time() - t
    1/0

    Genome.save_bed(query)
    1/0

    """
    for transcript in g.knownGene:
        print transcript, transcript.sequence()[:100] + "..."
        if transcript.txEnd > 8000: break
    """

    kg = g.table('knownGene')
    q = kg.select(kg.c.txStart < 5000)

    print list(g.session.execute(q))
