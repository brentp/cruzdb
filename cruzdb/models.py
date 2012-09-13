from sqlalchemy import Column, String
from sqlalchemy.ext.declarative import declared_attr
# needed to avoid circular imports
#CHANGED:from init import Base
from sequence import sequence as _sequence
from __init__ import Genome

"""
This is used to create a Model with the appropriate methods
from a UCSC table. It uses sqlalchemy reflection to
do the lifiting.

"""

class CruzException(Exception):
    pass




class ABase(object):
    _prefix_chain = ("tx", "chrom")
    @declared_attr
    def __tablename__(cls):
        return cls.__name__
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
        if not self.is_gene_pred: return []
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
        if len(cds) == 0: return []
        # grab all the sequences at once to reduce number of requests.
        all_seq = _sequence(self.db, self.chrom, cds[0][0] + 1, cds[-1][1])
        if len(cds) == 1:
            return all_seq
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
        """
        return the bins for efficient querying
        """
        return Genome.bins(self.start, self.end)

    def _introns(self, exons=None):
        if not self.is_gene_pred: return []
        se = self.exons
        if not (se or exons) or exons == []: return []
        starts, ends = zip(*exons) if exons is not None else zip(*se)
        return [(e, s) for e, s in zip(ends[:-1], starts[1:])]

    introns = property(_introns)

    def _xstream(self, s, e):
        f = Feature()
        f.txStart = s
        f.txEnd = e
        f.name = "region"
        f.cdsStart = f.cdsEnd = s
        f.strand = self.strand
        f.chrom = self.chrom
        return f

    def is_upstream_of(self, other):
        if self.chrom != other.chrom: return None
        if getattr(other, "strand", None) != "-":
            return self.end <= other.start
        # other feature is on - strand, so this must have higher start
        return self.start >= other.end

    def is_downstream_of(self, other):
        if self.chrom != other.chrom: return None
        if getattr(other, "strand", None) != "-":
            return self.start >= other.end
        # other feature is on - strand, so this must have higher start
        return self.end <= other.start

    def upstream(self, distance):
        """
        return the (start, end) of the region before the geneStart
        """
        if getattr(self, "strand", None) == "+":
            e = self.start
            s = e - distance
        else:
            s = self.end
            e = s + distance
        return self._xstream(s, e)

    def downstream(self, distance):
        """
        return the (start, end) of the region before the geneStart
        """
        if getattr(self, "strand", None) == "+":
            s = self.end
            e = s + distance
        else:
            e = self.start
            s = e - distance
        return self._xstream(s, e)

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

    def bed(self, *attrs, **kwargs):
        if self.is_gene_pred:
            return self.bed12(**kwargs)
        return "\t".join(map(str, (
                             [self.chrom, self.start, self.end] +
                             [getattr(self, attr) for attr in attrs]
                         )))



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
        name = self.name2 + "," + self.name if hasattr(self, "name2") \
                                            else self.name
        return "\t".join(map(str, (
            self.chrom, self.txStart, self.txEnd, name,
            score, self.strand, self.cdsStart, self.cdsEnd, rgb,
            len(exons), sizes, starts)))

    def localize(self, *positions, **kwargs):
        """
        convert global coordinate(s) to local taking
        introns into account and cds/tx-Start deinding on cdna=True kwarg
        """
        cdna = kwargs.get('cdna', False)
        # TODO: account for strand ?? add kwarg ??
        # if it's to the CDNA, then it's based on the cdsStart
        start, end = (self.cdsStart, self.cdsEnd) if cdna else \
                                        (self.start, self.end)
        introns = self.introns or None
        if cdna:
            if not self.is_coding:
                return ([None] * len(positions)) if len(positions) > 1 else None
            introns = self._introns(self.cds) or None

        if introns is None:
            print start, positions, end
            local_ps = [p - start if (start <= p < end) else None for p in positions]
            return local_ps[0] if len(positions) == 1 else local_ps

        introns = [(s - start, e - start) for s, e in introns]
        positions = [p - start for p in positions]
        # now both introns and positions are local starts based on cds/tx-Start
        local_ps = []
        l = end - start
        for original_p in positions:
            subtract = 0
            p = original_p
            if p < 0 or p >= l: # outside of transcript
                local_ps.append(None)
                continue
            for s, e in introns:
                # within intron
                if s <= p <= e:
                    subtract = None
                    break
                # otherwise, adjust for intron length.
                elif p >= e:
                    subtract += (e - s)

            local_ps.append(p - subtract if subtract is not None else None)

        assert all(p >=0 or p is None for p in local_ps), (local_ps)
        return local_ps[0] if len(positions) == 1 else local_ps

class Mixin(ABase):
    name = Column(String, unique=True, primary_key=True)

class kgXref(ABase):
    kgID = Column(String, unique=True, primary_key=True)

    def __repr__(self):
        return "%s(%s)" % (self.__tablename__, self.kgID)



Feature = Mixin
