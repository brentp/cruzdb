import unittest
from cruzdb import Genome
import os


class TestFeature(unittest.TestCase):
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

    def test_blat(self):
        try:
            import requests
        except ImportError:
            return
        g = Genome('hg18')
        f = g.refGene[19]
        f.chrom = "chr6"
        f.txStart = 135646802
        f.txEnd = 135646832
        r = list(f.blat())
        self.assert_(str(f.txStart) in repr(r), r)
        self.assert_(str(f.txEnd) in repr(r), r)

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
        query = g.knownGene.filter(and_(g.knownGene.txStart > 10000, g.knownGene.txEnd < 20000))
        c = StringIO()
        Genome.save_bed(query, c)
        c.seek(0)
        rows = c.readlines()
        for toks in (row.split("\t") for row in rows):
            self.assert_(len(toks) == 12)
            self.assert_(int(toks[1]) > 10000)
            self.assert_(int(toks[2]) < 20000)

    def test_link(self):
        feat = self.db.knownGene.first()
        l = feat.link()

        self.assert_(l == "http://genome.ucsc.edu/cgi-bin/hgGene?hgg_gene=uc001aaa.2&db=hg18", l)



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

class TestDb(unittest.TestCase):
    def setUp(self):
        self.dba = Genome('hg18', host="localhost", user="brentp")
        self.dbb = Genome('hg19', host="localhost", user="brentp")

    def test_protein(self):

        g = self.dba.knownGene.filter_by(name="uc010ntk.1").first()
        prot = g.protein

        self.assert_(prot.startswith("MIITQTSHCYMTSLGILFLINILPGTTGQGESRRQEPGDFVKQDIG"), prot)
        self.assert_(prot.endswith("SAIKGMIRKQ"), prot)

    def test_ok(self):
        ga = self.dba.refGene.filter_by(name2="MUC5B").first()
        self.assert_(ga is not None)

    def test_bins(self):
        bins = Genome.bins(12345, 56779)
        expected = set([1, 9, 73, 585])
        self.assertEqual(bins, expected)

    def test_tables(self):
        self.dba.refGene
        self.assert_("refGene" in self.dba.tables, self.dba.tables)

    def test_nearest(self):
        from cruzdb.models import Feature
        f = Feature()
        f.chrom = "chr1"
        f.txStart = 10
        f.txEnd = 61
        #db = Genome('hg18', host="localhost", user="brentp")
        db = self.dba
        self.assert_(db.refGene.first() is not None)
        self.assert_(db.refGene is not None)

        res = db.knearest(db.refGene, f, k=2)
        self.assert_(len(res) >= 2)

        f = db.refGene.first()
        key = (f.chrom, f.start, f.end, f.name)

        for k in (2, 4, 6):
            res = db.knearest("refGene", f, k=k)
            assert len(res) >= k
            self.assert_(key in ((n.chrom, n.start, n.end, n.name) for n in res),
                    (res, f))


        f = db.refGene.order_by(db.refGene.txStart).filter(db.refGene.c.strand == "+").first()
        assert f in db.upstream(db.refGene, f)

        down = db.downstream(db.refGene, f, k=10)
        self.assert_(len(down) >= 10)

        self.assert_(all(d.start >= f.start for d in down))


    def test_down_neg(self):
        db = self.dba
        fm = db.refGene.filter(db.refGene.c.strand == "-").first()
        down = db.downstream(db.refGene, fm, k=10)

        self.assert_(all(d.start <= fm.start for d in down))

    def test_dataframe(self):
        kg = Genome('hg18', host='localhost').dataframe('knownGene', limit=10)
        assert kg.shape[0] == 10

    def test_mirror(self):

        try:
            os.unlink('/tmp/__u.db')
        except OSError:
            pass
        g = Genome('hg18')
        g.mirror(['chromInfo'], 'sqlite:////tmp/__u.db')
        a = str(g.chromInfo.filter().first())

        gs = Genome('sqlite:////tmp/__u.db')

        b = str(gs.chromInfo.filter().first())
        self.assertEqual(a, b)
        os.unlink('/tmp/__u.db')

    def tearDown(self):
        del self.dba

if __name__ == "__main__":
    unittest.main()
