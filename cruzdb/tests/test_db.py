import unittest
from cruzdb import Genome


class TestMixin(unittest.TestCase):
    """
    class to just use the mixin without a connection
    """
    def setUp(self):
        from cruzdb.models import Feature
        self.f = Feature()
        self.f.chrom = "chr1"
        self.f.txStart = 10
        self.f.txEnd = 61

        self.f.cdsStart = 29
        self.f.cdsEnd =  59
        """
        + exon
        | coding-exon
        _ UTR
        - intron

        10        20    26 29   34   39      47   52     59 61
        ++++++++++______+++|||||-----||||||||-----|||||||+++

        # introns = [(20, 26), (34, 39), (47, 52)]

        # coding introns = [(34, 39), (47, 52)]
        """
        self.f.exonStarts = "10,26,39,52,"
        self.f.exonEnds = "20,34,47,61,"

        self.strand = self.f.strand = '+'


    def test_localize_out_of_bounds(self):
        f = self.f
        #self.assertEqual(f.localize(0, 80), [None, None])
        self.assertEqual(f.localize(0, 80, 61, 60, cdna=True),
                                     [None, None, None, None])

        self.assertEqual(f.localize(0, 80, 61, 60, cdna=False),
                                    [None, None, None, 34])

    def test_localize_in_intron(self):
        f = self.f
        self.assertEqual(f.localize(34, cdna=True), None)

    def test_localize_cdsBounds(self):
        f = self.f
        self.assertEqual(f.localize(f.cdsStart, cdna=True), 0)
        # DOES THIS MAKE SENSE? if it's at cdsEnd, it's None
        self.assertEqual(f.localize(f.cdsEnd, cdna=True), None)

        l = sum(e - s for s, e in f.cds)
        self.assertEqual(f.localize(f.cdsEnd - 1, cdna=True), l - 1)

    def test_localize_txBounds(self):
        f = self.f
        self.assertEqual(f.localize(f.txStart, cdna=False), 0)
        # DOES THIS MAKE SENSE? if it's at cdsEnd, it's None
        self.assertEqual(f.localize(f.txEnd, cdna=False), None)

        l = sum(e - s for s, e in f.exons)
        self.assertEqual(f.localize(f.txEnd - 1, cdna=False), l - 1)

    def test_upstream(self):
        f = self.f
        u = f.upstream(10)
        self.assertEqual(u.end, f.start)
        self.assert_(u.is_upstream_of(f))

    def test_downstream(self):
        f = self.f
        u = f.downstream(10)
        self.assertEqual(f.end, u.start)
        self.assert_(u.is_downstream_of(f))
        u.chrom = "fake"
        self.assertEqual(None,  u.is_upstream_of(f))

class TestBasic(unittest.TestCase):
    def setUp(self):
        self.db = Genome('hg18', host="localhost", user="brentp")

    def testFirst(self):
        self.assert_(hasattr(self.db.refGene.first(), "txStart"))

    def test_bed_gene_pred(self):
        g = Genome('hg19', host="localhost", user="brentp")
        from sqlalchemy import and_
        from cStringIO import StringIO
        query = g.knownGene.filter(and_(g.table('knownGene').c.txStart > 10000, g.table('knownGene').c.txEnd < 20000))
        c = StringIO()
        Genome.save_bed(query, c)
        c.seek(0)
        rows = c.readlines()
        for toks in (row.split("\t") for row in rows):
            self.assert_(len(toks) == 12)
            self.assert_(int(toks[1]) > 10000)
            self.assert_(int(toks[2]) < 20000)

    def test_bed_other(self):
        g = self.db
        self.assertEqual(g.cpgIslandExt[12].bed(), 'chr1	829557	830482')
        self.assertEqual(g.cpgIslandExt[12].bed('length', 'perCpg'), 'chr1	829557	830482	925	17.9')


class TestGene(unittest.TestCase):
    def setUp(self):
        self.db = Genome('hg18', host="localhost", user="brentp")
        self.gene = self.db.refGene.filter_by(name2="MUC5B").first()

    def testExons(self):
        self.assert_(isinstance(self.gene.exons, list))
        self.assert_(self.gene.exons[0][0] >= self.gene.txStart)
        self.assert_(self.gene.exons[0][0] <= self.gene.cdsStart)

    def testIntrons(self):
        self.assert_(isinstance(self.gene.introns, list))
        self.assert_(self.gene.introns[0][0] >= self.gene.txStart)
        self.assert_(self.gene.introns[0][0] == self.gene.exons[0][1])
        self.assert_(all((s < e) for s, e in self.gene.introns))

    def testBed12(self):
        expected = "chr1	891739	900347	PLEKHN1,NM_032129	0	+	891774	899818	.	16	118,100,147,81,73,128,96,81,76,137,150,141,141,219,49,663,	0,207,3780,4024,4189,4382,4616,4827,5578,5791,6364,6689,7003,7336,7819,7945,"
        transcript = self.db.refGene.filter_by(name="NM_032129").first()
        self.assertEqual(transcript.bed12(), expected)

class Test2Db(unittest.TestCase):
    def setUp(self):
        self.dba = Genome('hg18', host="localhost", user="brentp")
        self.dbb = Genome('hg19', host="localhost", user="brentp")

    def test_ok(self):
        ga = self.dba.refGene.filter_by(name2="MUC5B").first()
        gb = self.dbb.refGene.filter_by(name2="MUC5B").first()
        assert ga.txStart != gb.txStart
        # TODO: these shouldb e different, but they both hg18
        # session stuff?
        print ga.db, gb.db

    def test_bins(self):
        bins = Genome.bins(12345, 56779)
        expected = set([1, 9, 73, 585])
        self.assertEqual(bins, expected)

if __name__ == "__main__":
    unittest.main()
