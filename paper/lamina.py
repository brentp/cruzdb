import os.path as op
from toolshed import reader
from cruzdb import Genome

def lamina():
    if not op.exists('lamina.bed'):
        fh = open('lamina.bed', 'w')
        fh.write("#chrom\tstart\tend\tvalue\n")
        for gff in reader('http://www.nature.com/nature/journal/v453/n7197/extref/nature06947-s2.txt', header=False):
            fh.write("\t".join([gff[0], gff[3], gff[4], gff[5]]) + "\n")
        fh.close()
    return 'lamina.bed'

fname = 'supplement/Additional-File-11_lamina.anno.bed'
hg18 = Genome('sqlite:///hg18.db')
if not op.exists(fname):
    fhout = open(fname, 'w')
    hg18.annotate(lamina(), ('refGene', ), feature_strand=True, in_memory=True, parallel=True, out=fhout)
    fhout.close()


for cutoff in (0.90, 0.95):
    fh = open('/tmp/genes-%.2f.txt' % cutoff, 'w')
    for d in reader(fname):
        if float(d['value']) < cutoff: continue
        if d['refGene_distance'] == '0' or \
           d['refGene_distance'].startswith("0;"):
            print >>fh, "\n".join(d['refGene_name'].split(";"))
    fh.close()

cutoff = 0.90
fh = open('/tmp/genes-overlap-complete.txt', 'w')
for d in (l for l in reader(lamina()) if float(l['value']) > cutoff):
    if float(d['value']) < cutoff: continue

    start, end = map(int, (d['start'], d['end']))

    res = hg18.bin_query('refGene', d['chrom'], start, end).all()

    if len(res) == 0: continue

    for r in res:
        # genes completely contained within an LAD
        if start <= r.start and end >= r.end:
            print >>fh, r.gene_name
