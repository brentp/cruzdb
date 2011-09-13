import unittest
from cruzdb import Genome


class TestBasic(unittest.TestCase):
    def setUp(self):
        self.db = Genome('hg18')

    def testFirst(self):
        self.assert_(hasattr(self.db.refGene.first(), "txStart"))

    def test_bed_gene_pred(self):
        g = Genome('hg19')
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
        self.db = Genome('hg18')
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
        self.dba = Genome('hg18')
        self.dbb = Genome('hg19')

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
