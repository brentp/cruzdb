from init import initialize_sql
from sqlalchemy import create_engine
import sys

from tests import test

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
        self.models = __import__("models", globals(), locals(), ["Mixin"], -1)

    def _map(self, table):
        # if the table hasn't been mapped, do so here.
        #if not table in self.models.Base.metadata.tables:
        if not table in self.__tables:
            # make a new class
            self.__tables[table] = type(table, (self.Base, self.models.Mixin), {})

    def __getattr__(self, table):
        self._map(table)
        return self.session.query(self.__tables[table])

    def table(self, table):
        self._map(table)
        return self.Base.metadata.tables[table]

    def sql(self, query):
        return self.engine.execute(query)

    def __repr__(self):
        return "%s('%s')" % (self.__class__.__name__,
            self.url % dict(db=self.db, user=self.user, host=self.host,
            password="?" if self.password else ""))

    def bed(self, filename=sys.stdout):
        pass
        # TODO allow writing bed6 and bed12...

    def mirror(self, tables, to_url):
        to_engine = create_engine(to_url)


if __name__ == "__main__":
    g = Genome(db="hg18", host="localhost", user="")


    import sys
    sys.exit()
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

    print g.knownGene.filter(g.table('knownGene').c.txStart > 10000).first()

    """
    for transcript in g.knownGene:
        print transcript, transcript.sequence()[:100] + "..."
        if transcript.txEnd > 8000: break
    """

    kg = g.table('knownGene')
    q = kg.select(kg.c.txStart < 5000)

    print list(g.session.execute(q))

    print g.cpgIslandExt[12]



#################################

