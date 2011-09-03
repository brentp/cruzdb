from sqlalchemy import Column, String
from sqlalchemy.ext.declarative import declared_attr
# needed to avoid circular imports
#CHANGED:from init import Base
from sequence import sequence as _sequence

"""
This is used to create a Model with the appropriate methods
from a UCSC table. It uses sqlalchemy reflection to
do the lifiting.

"""

class CruzException(Exception):
    pass


class Mixin(object):
    _prefix_chain = ("tx", "chrom")

    @declared_attr
    def __tablename__(cls):
        return cls.__name__
    name = Column(String, unique=True, primary_key=True)
    __table_args__ = {'autoload': True}

    @property
    def is_coding(self):
        try:
            return self.cdsStart != self.cdsEnd
        except AttributeError:
            return False

    @property
    def exons(self):
        # drop the trailing comma
        starts = (int(s) for s in self.exonStarts[:-1].split(","))
        ends = (int(s) for s in self.exonEnds[:-1].split(","))
        return zip(starts, ends)

    @property
    def coding_exons(self):
        """
        includes the entire exon as long as any of it is > cdsStart and <
        cdsEnd
        """
        # drop the trailing comma
        starts = (int(s) for s in self.exonStarts[:-1].split(","))
        ends = (int(s) for s in self.exonEnds[:-1].split(","))
        return [(s, e) for s, e in zip(starts, ends)
                                          if e > self.cdsStart and
                                             s < self.cdsEnd]

    @property
    def cds(self):
        """just the parts of the exons that are translated"""
        ces = self.coding_exons
        if len(ces) < 1: return ces
        ces[0] = (self.cdsStart, ces[0][1])
        ces[-1] = (ces[-1][0], self.cdsEnd)
        assert all((s < e for s, e in ces))
        return ces

    @property
    def cds_sequence(self):
        seqs = []
        cds = self.cds
        # grab all the sequences at once to reduce number of requests.
        all_seq = _sequence(self.db, self.chrom, cds[0][0] + 1, cds[-1][1])
        lowest = cds[0][0]
        cds0 = [(s - lowest, e - lowest) for s, e in cds]
        for cstart, cend in cds0:
            seqs.append(all_seq[cstart:cend])
        return seqs

    @property
    def browser_link(self):
        return "http://genome.ucsc.edu/cgi-bin/hgTracks?db=%s&position=%s" % (self.db, self.position)

    @property
    def position(self):
        return "%s:%i-%i" % (self.chrom, self.start, self.end)

    @property
    def bins(self):
        bins = [1];
        for b in range(1 + (start>>26), 1 + ((end-1)>>26)+1):     bins.append(b)
        for b in range(9 + (start>>23), 9 + ((end-1)>>23)+1):     bins.append(b)
        for b in range(73 + (start>>20), 73 + ((end-1)>>20)+1):   bins.append(b)
        for b in range(585 + (start>>17), 585 + ((end-1)>>17)+1): bins.append(b)
        return set(bins)

    @property
    def introns(self):
        starts, ends = zip(*self.exons)
        return [(e, s) for e, s in zip(ends[:-1], starts[1:])]

    @property
    def utr5(self):
        if not self.is_coding or len(self.exons) < 2: return (None, None)
        if self.strand == "+":
            return (self.txStart, self.cdsStart)
        else:
            return (self.cdsEnd, self.txEnd)

    @property
    def utr3(self):
        if not self.is_coding or len(self.exons) < 2: return (None, None)
        if self.strand == "-":
            return (self.txStart, self.cdsStart)
        else:
            return (self.cdsEnd, self.txEnd)

    def __len__(self):
        return self.end - self.start

    def __cmp__(self, other):
        if self.chrom != other.chrom: return 0
        if self.start < other.start: return -1
        return 1

    @property
    def start(self):
        for prefix in self._prefix_chain:
            try: return getattr(self, prefix + "Start")
            except AttributeError: pass
        raise Exception("no start for %r" % self)

    @property
    def end(self):
        for prefix in self._prefix_chain:
            try: return getattr(self, prefix + "End")
            except AttributeError: pass
        raise Exception("no end for %r" % self)


    def __repr__(self):
        return "%s(%s:%s:%i-%i)" % (self.__tablename__, self.chrom, self.gene_name,
                self.start, self.end)

    @property
    def gene_name(self):
        if hasattr(self, "name2"): return self.name2
        if hasattr(self, "name"): return self.name
        return self.position

    @property
    def db(self):
        # grab the database name from the current row
        # e.g. hg18
        return self.metadata.bind.url.database

    def __str__(self):
        # output something bed-like
        # TODO: for most tables, this could be bed12
        fields = "chrom start end gene_name".split()
        s = "\t".join(map(str, (getattr(self, field) for field in fields)))
        if hasattr(self, "strand"):
            s += "\t.\t%s" % (self.strand)
        return s

    def sequence(self, per_exon=False):
        """
        if per-exon is True, return an array of exon sequences
        NOTE: this is never reverse-complemented. TODO??
        """
        db = self.db
        if not per_exon:
            start = self.txStart + 1
            return _sequence(db, self.chrom, start, self.txEnd)
        else:
            # TODO: use same strategy as cds_sequence to reduce # of requests.
            seqs = []
            for start, end in self.exons:
                seqs.append(_sequence(db, self.chrom, start + 1, end))
            return seqs

    @property
    def is_gene_pred(self):
        """
        http://genome.ucsc.edu/FAQ/FAQformat.html#format9
        """
        return hasattr(self, "exonStarts")

    def bed12(self, score="0", rgb="."):
        """convert the exon stuff into bed12
        http://genome.ucsc.edu/FAQ/FAQformat.html#format1
        # TODO: all name_fn kwarg like:
            name_fn = lambda feat: feat.name + "," feat.name2
        in case we want to show gene name and transcript name.
        """
        if not self.is_gene_pred:
            raise CruzException("can't create bed12 from non genepred feature")
        exons = self.exons
        # go from global start, stop, to relative start, length...
        sizes = ",".join([str(e[1] - e[0]) for e in exons]) + ","
        starts = ",".join([str(e[0] - self.txStart) for e in exons]) + ","
        return "\t".join(map(str, (
            self.chrom, self.txStart, self.txEnd, self.name2 + "," + self.name,
            score, self.strand, self.cdsStart, self.cdsEnd, rgb,
            len(exons), sizes, starts)))




