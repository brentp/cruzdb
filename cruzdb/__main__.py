
def main(input_bed, genome_version, tables, feature_strand=True):
    from . import Genome
    Genome(genome_version).annotate(input_bed, tables,
            feature_strand=feature_strand)

if __name__ == "__main__":
    import sys
    tables = sys.argv[3:] if len(sys.argv) > 3 else ('refGene', 'cpgIslandExt')
    main(sys.argv[1], sys.argv[2], tables)
