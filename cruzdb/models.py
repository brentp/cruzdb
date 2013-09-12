"""
This is used to create a Model with the appropriate methods
from a UCSC table. It uses sqlalchemy reflection to
do the lifiting.

"""
from sqlalchemy import Column, String, ForeignKey, Float, Integer
from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy.orm import relationship, backref
from sqlalchemy.schema import PrimaryKeyConstraint

import sys
from operator import itemgetter

# needed to avoid circular imports
#CHANGED:from init import Base
from sequence import sequence as _sequence
from __init__ import Genome


import re

def _ncbi_parse(html):
    from collections import OrderedDict

    try:
        info = html.split("Sequences producing significant alignments")[1].split("<tbody>")[1]
    except IndexError:
        print >>sys.stderr, html
        raise
    info = info.split("</table>")[0]
    regexp = re.compile(r'<tr>(.+?)(<\/tr>)', re.MULTILINE | re.DOTALL)
    tdreg = re.compile(r'<td.*?>(.+?)(?:</td>)', re.MULTILINE | re.DOTALL)
    colnames = ("accession", "org", "description", "max_score", "total_score",
                    "query_coverage", "e_value", "max_ident", "link")
    for record in (r.groups(0)[0] for r in regexp.finditer(info)):
        try:
            cols = tdreg.findall(record)
            pcols = [c.split(">")[1].split("<")[0].strip() if "<" in c else c.strip() for c in cols[:-1]]
            pcols.insert(1, " ".join(pcols[1].replace("PREDICTED: ", "").split(" ")[:2]))
            try:
                pcols.append(cols[-1].split("href=")[1].split(">")[0])
            except IndexError: # no link
                pcols.append("")
            yield OrderedDict(zip(colnames, pcols))
        except:
            print >>sys.stderr, record

class CruzException(Exception):
    pass

class Interval(object):
    """
    Interval class for convenience

    Parameters
    ----------

    start : int
    
    end : int

    chrom : str

    name : str
        optional name for the interval
    """

    __slots__ = ('chrom', 'start', 'end', 'name', 'gene_name')
    def __init__(self, start, end, chrom=None, name=None):
        self.start, self.end = start, end
        self.chrom = chrom
        self.name = self.gene_name = name

    def overlaps(self, other):
        """
        check for overlap with the other interval
        """
        if self.chrom != other.chrom: return False
        if self.start > other.end: return False
        if other.start > self.end: return False
        return True

    def is_upstream_of(self, other):
        """
        check if this is upstream of the `other` interval taking the strand of
        the other interval into account
        """
        if self.chrom != other.chrom: return None
        if getattr(other, "strand", None) == "+":
            return self.end <= other.start
        # other feature is on - strand, so this must have higher start
        return self.start >= other.end

    def distance(self, other_or_start=None, end=None, features=False):
        """
        check the distance between this an another interval
        Parameters
        ----------

        other_or_start : Interval or int
            either an integer or an Interval with a start attribute indicating
            the start of the interval

        end : int
            if `other_or_start` is an integer, this must be an integer
            indicating the end of the interval

        features : bool
            if True, the features, such as CDS, intron, etc. that this feature
            overlaps are returned.
        """
        if end is None:
            assert other_or_start.chrom == self.chrom

        other_start, other_end = get_start_end(other_or_start, end)

        if other_start > self.end:
            return other_start - self.end
        if self.start > other_end:
            return self.start - other_end
        return 0

class ABase(object):
    """
    Base object that wraps returned database rows
    """

    _prefix_chain = ("tx", "chrom")
    @declared_attr
    def __tablename__(cls):
        return cls.__name__
    __table_args__ = {'autoload': True}
    __mapper_args__= {'always_refresh': False, 'exclude_properties': ['dist', '_dist']}

    anno_cols = ("name", "distance", "feature")

    @property
    def is_coding(self):
        try:
            return self.cdsStart != self.cdsEnd
        except AttributeError:
            return False

    def _repr_html_(self):
        info = dict(db=self.db, start=self.start, end=self.end,
                chrom=self.chrom)

        info['name'] = self.__class__.name
        url = "http://genome.ucsc.edu/cgi-bin/hgTracks?position=%(chrom)s:%(start)s-%(end)s&db=%(db)s&%(name)s=full"  % info
        return "<iframe src='%s' style='margin:0px; border:0; height:500px; width:100%%' ></iframe>" % url

    @property
    def exons(self):
        """
        return a list of exons [(start, stop)] for this object if appropriate
        """
        # drop the trailing comma
        if not self.is_gene_pred: return []
        if hasattr(self, "exonStarts"):
            starts = (long(s) for s in self.exonStarts[:-1].split(","))
            ends = (long(s) for s in self.exonEnds[:-1].split(","))
        else: # it is bed12
            starts = [self.start + long(s) for s in self.chromStarts[:-1].split(",")]
            ends = [starts[i] + long(size) for i, size \
                    in enumerate(self.blockSizes[:-1].split(","))]


        return zip(starts, ends)

    @property
    def gene_features(self):
        """
        return a list of features for the gene features of this object.
        This would include exons, introns, utrs, etc.
        """
        nm, strand = self.gene_name, self.strand
        feats = [(self.chrom, self.start, self.end, nm, strand, 'gene')]
        for feat in ('introns', 'exons', 'utr5', 'utr3', 'cdss'):
            fname = feat[:-1] if feat[-1] == 's' else feat
            res = getattr(self, feat)
            if res is None or all(r is None for r in res): continue
            if not isinstance(res, list): res = [res]
            feats.extend((self.chrom, s, e, nm, strand, fname) for s, e in res)

        tss = self.tss(down=1)
        if tss is not None:
            feats.append((self.chrom, tss[0], tss[1], nm, strand, 'tss'))
            prom = self.promoter()
            feats.append((self.chrom, prom[0], prom[1], nm, strand, 'promoter'))

        return sorted(feats, key=itemgetter(1))

    def tss(self, up=0, down=0):
        """
        Return a start, end tuple of positions around the transcription-start
        site

        Parameters
        ----------

        up : int
           if greature than 0, the strand is used to add this many upstream
           bases in the appropriate direction

        down : int
           if greature than 0, the strand is used to add this many downstream
           bases into the gene.
        """
        if not self.is_gene_pred: return None
        tss = self.txEnd if self.strand == '-' else self.txStart
        start, end = tss, tss
        if self.strand == '+':
            start -= up
            end += down
        else:
            start += up
            end -= down
            start, end = end, start
        return max(0, start), max(end, start, 0)

    def promoter(self, up=2000, down=0):
        """
        Return a start, end tuple of positions for the promoter region of this
        gene

        Parameters
        ----------

        up : int
           this distance upstream that is considered the promoter

        down : int
           the strand is used to add this many downstream bases into the gene.
        """
        if not self.is_gene_pred: return None
        return self.tss(up=up, down=down)

    @property
    def coding_exons(self):
        """
        includes the entire exon as long as any of it is > cdsStart and <
        cdsEnd
        """
        # drop the trailing comma
        starts = (long(s) for s in self.exonStarts[:-1].split(","))
        ends = (long(s) for s in self.exonEnds[:-1].split(","))
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

    cdss = cds

    def _cds_sequence(self, cds):
        seqs = []
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
    def cds_sequence(self):
        """
        a list of genomic sequences for the CDS's
        """
        return self._cds_sequence(self.cds)

    @property
    def mrna_sequence(self):
        """
        a list of genomic sequences for the mRNA's
        """
        return self._cds_sequence(self.coding_exons)

    @property
    def browser_link(self):
        return "http://genome.ucsc.edu/cgi-bin/hgTracks?db=%s&position=%s" % (self.db, self.position)

    @property
    def position(self):
        " a chrom:start-stop representation of this feature"
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
        """
        return a boolean indicating whether this feature is upstream of `other`
        taking the strand of other into account
        """
        if self.chrom != other.chrom: return None
        if getattr(other, "strand", None) == "-":
            # other feature is on - strand, so this must have higher start
            return self.start >= other.end
        return self.end <= other.start

    def is_downstream_of(self, other):
        """
        return a boolean indicating whether this feature is downstream of
        `other` taking the strand of other into account
        """
        if self.chrom != other.chrom: return None
        if getattr(other, "strand", None) == "-":
            # other feature is on - strand, so this must have higher start
            return self.end <= other.start
        return self.start >= other.end

    def features(self, other_start, other_end):
        """
        return e.g. "intron;exon" if the other_start, end overlap introns and
        exons
        """
        # completely encases gene.
        if other_start <= self.start and other_end >= self.end:
            return ['gene' if self.cdsStart != self.cdsEnd else 'nc_gene']
        other = Interval(other_start, other_end)
        ovls = []
        tx = 'txEnd' if self.strand == "-" else 'txStart'
        if hasattr(self, tx) and other_start <= getattr(self, tx) <= other_end \
            and self.cdsStart != self.cdsEnd:
                ovls = ["TSS"]
        for ftype in ('introns', 'exons', 'utr5', 'utr3', 'cdss'):
            feats = getattr(self, ftype)
            if not isinstance(feats, list): feats = [feats]
            if any(Interval(f[0], f[1]).overlaps(other) for f in feats):
                ovls.append(ftype[:-1] if ftype[-1] == 's' else ftype)
        if 'cds' in ovls:
            ovls = [ft for ft in ovls if ft != 'exon']
        if self.cdsStart == self.cdsEnd:
            ovls = ['nc_' + ft for ft in ovls]
        return ovls

    def distance(self, other_or_start=None, end=None, features=False):
        """
        check the distance between this an another interval
        Parameters
        ----------

        other_or_start : Interval or int
            either an integer or an Interval with a start attribute indicating
            the start of the interval

        end : int
            if `other_or_start` is an integer, this must be an integer
            indicating the end of the interval

        features : bool
            if True, the features, such as CDS, intron, etc. that this feature
            overlaps are returned.
        """
        if end is None:
            assert other_or_start.chrom == self.chrom

        other_start, other_end = get_start_end(other_or_start, end)

        if other_start > self.end:
            return other_start - self.end, "intergenic"
        if self.start > other_end:
            return self.start - other_end, "intergenic"
        if features: return (0, "+".join(self.features(other_start, other_end)))
        return (0, "")

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
        """
        return the 5' UTR if appropriate
        """
        if not self.is_coding or len(self.exons) < 2: return (None, None)
        if self.strand == "+":
            s, e = (self.txStart, self.cdsStart)
        else:
            s, e = (self.cdsEnd, self.txEnd)
        if s == e: return (None, None)
        return s, e

    @property
    def utr3(self):
        """
        return the 3' UTR if appropriate
        """
        if not self.is_coding or len(self.exons) < 2: return (None, None)
        if self.strand == "-":
            s, e = (self.txStart, self.cdsStart)
        else:
            s, e = (self.cdsEnd, self.txEnd)
        if s == e: return (None, None)
        return s, e

    def __len__(self):
        return self.end - self.start

    def __cmp__(self, other):
        if self.chrom != getattr(other, "chrom", other): return 0
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
        try:
            self.start
            return "%s(%s:%s:%i-%i)" % (self.__class__.__name__, self.chrom, self.gene_name,
                self.start, self.end)
        except:
            return "%s(%s)" % (self.__tablename__, self.chrom)


    @property
    def gene_name(self):
        if hasattr(self, "name2"): return self.name2
        if hasattr(self, "name"): return self.name
        return self.position

    @property
    def db(self):
        # grab the database name from the current row
        # e.g. hg18
        return self._table.bind.url.database

    def __str__(self):
        # output something bed-like
        fields = "chrom start end gene_name".split()
        s = "\t".join(map(str, (getattr(self, field) for field in fields)))
        if hasattr(self, "score"):
            s += "\t%.2f" % (self.score)
            if hasattr(self, "strand"):
                s += "\t%s" % (self.strand)

        elif hasattr(self, "strand"):
            s += "\t.\t%s" % (self.strand)
        return s

    def sequence(self, per_exon=False):
        """
        Return the sequence for this feature.
        if per-exon is True, return an array of exon sequences
        This sequence is never reverse complemented
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

    def __iter__(self):
        for k in ('chrom', 'start', 'end', 'name', 'score', 'strand'):
            yield str(getattr(self, k, ""))

    def ncbi_blast(self, db="nr", megablast=True, sequence=None):
        """
        perform an NCBI blast against the sequence of this feature
        """
        import requests
        requests.defaults.max_retries = 4
        assert sequence in (None, "cds", "mrna")
        seq = self.sequence() if sequence is None else ("".join(self.cds_sequence if sequence == "cds" else self.mrna_sequence))
        r = requests.post('http://blast.ncbi.nlm.nih.gov/Blast.cgi',
                        timeout=20,
                        data=dict(
                            PROGRAM="blastn",
                            #EXPECT=2,
                            DESCRIPTIONS=100,
                            ALIGNMENTS=0,
                            FILTER="L", # low complexity
                            CMD="Put",
                            MEGABLAST=True,
                            DATABASE=db,
                            QUERY=">%s\n%s" % (self.name, seq)
                        )
                    )

        if not ("RID =" in r.text and "RTOE" in r.text):
            print >>sys.stderr, "no results"
            raise StopIteration
        rid = r.text.split("RID = ")[1].split("\n")[0]

        import time
        time.sleep(4)
        print >>sys.stderr, "checking..."
        r = requests.post('http://blast.ncbi.nlm.nih.gov/Blast.cgi',
                data=dict(RID=rid, format="Text",
                    DESCRIPTIONS=100,
                    DATABASE=db,
                    CMD="Get", ))
        while "Status=WAITING" in r.text:
            print >>sys.stderr, "checking..."
            time.sleep(10)
            r = requests.post('http://blast.ncbi.nlm.nih.gov/Blast.cgi',
                data=dict(RID=rid, format="Text",
                    CMD="Get", ))
        for rec in _ncbi_parse(r.text):
            yield rec

    def blat(self, db=None, sequence=None, seq_type="DNA"):
        """
        make a request to the genome-browsers BLAT interface
        sequence is one of None, "mrna", "cds"
        returns a list of features that are hits to this sequence.
        """
        from . blat_blast import blat, blat_all
        assert sequence in (None, "cds", "mrna")
        seq = self.sequence() if sequence is None else ("".join(self.cds_sequence if sequence == "cds" else self.mrna_sequence))
        if isinstance(db, (tuple, list)):
            return blat_all(seq, self.gene_name, db, seq_type)
        else:
            return blat(seq, self.gene_name, db or self.db, seq_type)

    @property
    def is_gene_pred(self):
        """
        http://genome.ucsc.edu/FAQ/FAQformat.html#format9
        """
        return hasattr(self, "exonStarts") or hasattr(self, 'chromStarts')

    def bed(self, *attrs, **kwargs):
        """
        return a bed formatted string of this feature
        """
        exclude = ("chrom", "start", "end", "txStart", "txEnd", "chromStart",
                "chromEnd")
        if self.is_gene_pred:
            return self.bed12(**kwargs)
        return "\t".join(map(str, (
                 [self.chrom, self.start, self.end] +
                 [getattr(self, attr) for attr in attrs if not attr in exclude]
                         )))


    def bed12(self, score="0", rgb="."):
        """
        return a bed12 (http://genome.ucsc.edu/FAQ/FAQformat.html#format1)
        representation of this interval
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

    def globalize(self, position, cdna=True):
        1/0
        start, end = (self.cdsStart, self.cdsEnd) if cdna else \
                                        (self.start, self.end)
        exons = self.exons or None
        pos = position + start
        if exons is None:
            return pos

        subtract = 0
        print >>sys.stderr, "exon lengths:", sum((ie - ib) for ib, ie in self.exons)
        for estart, eend in exons:
            if iend < pos:
                subtract += (iend - istart)
            elif istart < pos and iend > pos:
                subtract += (pos - istart)
            print >>sys.stderr, subtract, (istart, iend), pos
        return pos - subtract

    def localize(self, *positions, **kwargs):
        """
        convert global coordinate(s) to local taking
        introns into account and cds/tx-Start depending on cdna=True kwarg
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
            print >>sys.stderr, p, l
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


def get_start_end(other_or_start, end):
    if end is None:
        other_start, other_end = other_or_start.start, other_or_start.end
    else:
        other_start, other_end = other_or_start, end
    return other_start, other_end

class Feature(ABase):
    name = Column(String, unique=False, primary_key=True)

class cpgIslandExt(Feature):
    anno_cols = ("name", "distance", "feature")

    def distance(self, other_or_start=None, end=None, features="unused",
            shore_dist=3000):
        """
        check the distance between this an another interval
        Parameters
        ----------

        other_or_start : Interval or int
            either an integer or an Interval with a start attribute indicating
            the start of the interval

        end : int
            if `other_or_start` is an integer, this must be an integer
            indicating the end of the interval

        features : bool
            if True, the features, such as CDS, intron, etc. that this feature
            overlaps are returned.
        """
        # leave features kwarg to match signature from Feature.distance
        if end is None:
            assert other_or_start.chrom == self.chrom
        other_start, other_end = get_start_end(other_or_start, end)

        dist = 0
        if other_start > self.end:
            dist = other_start - self.end
        elif self.start > other_end:
            dist = self.start - other_end
        assert dist >= 0

        if dist > 0: dist = (dist, "shore" if abs(dist) <= shore_dist else "")
        else: dist = (0, "island")
        return dist

cpgRafaLab = cpgIslandExt


class SNP(ABase):
    __table_args__ = (
            PrimaryKeyConstraint('name', 'chrom', 'chromStart'),
            dict(autoload=True),)
    # can't add name or it gives error on select.
    #name = Column(String, unique=False)
    @property
    def name2(self):
        return self.name + (("-" + self.func) if self.func != "unknown" else "")

    def to_simple(self):
        return Interval(self.chromStart, self.chromEnd, self.chrom, self.name2)

class chromInfo(ABase):
    def __repr__(self):
        return "%s(%s:%i)" % (self.__tablename__, self.chrom, self.size)

    __str__ = __repr__



class Blat(Feature):

    identity = Column(Float)
    span = Column(Integer)
    db = Column(String(6))

    def __str__(self):
        res = Feature.__str__(self).replace("\t.\t", "\t%.1f%%\t" % self.identity)
        res += "\t%s\t%s" % (self.span, self.db)
        return res

    @property
    def score(self):
        return self.identity

    @property
    def hit_length(self):
        return self.span

class kgXref(ABase):
    __tablename__ = "kgXref"

    kgID = Column(String, primary_key=True)

#    @declared_attr
#    def kgID(self):
#        return Column(String, ForeignKey('knownGene.name'), primary_key=True)


    def __repr__(self):
        return "%s(%s/%s)" % (self.__tablename__, self.geneSymbol, self.kgID)

    __str__ = __repr__

class knownGene(ABase):
    __tablename__ = "knownGene"

    __mapper_args__= {'always_refresh': False, 'exclude_properties': ['dist',
        '_dist']}

    __preload_classes__ = ("kgXref",)

    anno_cols = ("name", "distance", "feature")

    @declared_attr
    def name(cls):
        return Column(ForeignKey('kgXref.kgID'), primary_key=True)

    #@declared_attr
    #def kgXref(cls):
    #    #return relationship("kgXref", backref=backref("knownGene"), lazy="subquery")
#
#        return relationship(lambda: kgXref,
#                primaryjoin=lambda: knownGene.name==kgXref.kgID,
#            lazy="subquery")
            #viewonly=True)

    @property
    def name2(self):
        return self.kgXref.geneSymbol

    def link(self):
        l = "http://genome.ucsc.edu/cgi-bin/hgGene?hgg_gene=%s&db=%s"
        return l % (self.name, self.db)

    @property
    def protein(self):
        from toolshed import nopen
        l = "http://genome.ucsc.edu/cgi-bin/hgGene?hgg_do_getProteinSeq=1&hgg_gene="
        url = l + self.name
        seq = [x.strip() for x in nopen(url) if x.strip() and
                not ">" in x]
        return "".join(seq)


