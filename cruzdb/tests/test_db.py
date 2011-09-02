import unittest
from cruzdb import Genome


class TestBasic(unittest.TestCase):
    def setUp(self):
        self.db = Genome('hg18')

    def testFirst(self):
        self.assert_(hasattr(self.db.refGene.first(), "txStart"))

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



if __name__ == "__main__":
    unittest.main()
