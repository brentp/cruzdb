from cruzdb import Genome
from cruzdb.sequence import sequence

# mirror the neede tables from UCSC to a local sqlite db
local = Genome('hg19').mirror(('refGene', 'targetScanS'), 'sqlite:///hg19.mirna.db')

# connect to the newly created local sqlite database instance.
refseq_ids = []

# iterate over the coding in refGene
for gene in (rgene for rgene in local.refGene if rgene.is_coding):

    if None in gene.utr3: continue # skip genes with no UTR

    utr_start, utr_end = gene.utr3
    # query the targetScan miRNA table with efficient bin query 
    sites = local.bin_query('targetScanS', gene.chrom, utr_start, utr_end)

    # print BED file of genes whose 3'UTR contains a miR-96 target site
    # with a score > 85.
    if any("miR-96" in s.name and s.score > 85 for s in sites):
        refseq_ids.append(gene.name) # save the refSeq for later GO analysis

        # gene is a python object but its string representation is BED format
        # we also print out the UTR sequence.
        print gene, sequence('hg19', gene.chrom, utr_start, utr_end)

# open a webbrowser to show enrichment of the genes we've selected in DAVID
Genome.david_go(refseq_ids)
