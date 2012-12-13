import requests
from .models import Blat

def blat(seq, name, db, seq_type="DNA"):
    r = requests.post('http://genome.ucsc.edu/cgi-bin/hgBlat',
            data=dict(db=db, type=seq_type, userSeq=seq, output="html"))
    if "Sorry, no matches found" in r.text:
        raise StopIteration
    text = r.text.split("<TT><PRE>")[1].split("</PRE></TT>")[0].strip().split("\n")
    for istart, line in enumerate(text):
        if "-----------" in line: break
    istart += 1
    for i, hit in enumerate(t.rstrip("\r\n") for t in text[istart:]):
        hit = hit.split(" YourSeq ")[1].split()
        f = Blat()
        # blat returns results without chr prefix
        if not hit[5].startswith("chr"): hit[5] = "chr" + hit[5]
        f.chrom = hit[5]
        f.txStart = long(hit[7]) - 1 # blat returns 1-based hits.
        f.txEnd = long(hit[8])
        f.strand = hit[6]
        f.identity = float(hit[4].rstrip("%"))
        f.span = int(hit[-1])
        f.db = db
        f.name = "blat-hit-%i-to-%s (%i bp)" % (i + 1, name, f.span)
        yield f

def blat_all(seq, name, dbs, seq_type="DNA"):
    for db in dbs:
        for f in blat(seq, name, db, seq_type):
            yield f


