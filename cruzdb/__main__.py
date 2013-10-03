"""
use cruzdb to annotate an input file with tables from the corresponding
database. example usage:

    $ python -m cruzdb hg19 intervals.bed refGene dgv snp137Common wgEncodeRegTfbsClusteredV2

To annotate with the 4 tables. Output is to stdout.
"""

def annotate(input_bed, genome_version, tables, feature_strand=True,
        in_memory=False):
    from . import Genome
    Genome(genome_version).annotate(input_bed, tables,
            feature_strand=feature_strand, in_memory=in_memory)

def main():
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("genome", help="genome version to use (e.g. hg19)")
    p.add_argument("bed", help="input bed file to annotated")
    p.add_argument("--in-memory", default=False, action="store_true")

    p.add_argument('tables', nargs='*', help='tables to annotated with',
            default=['refGene', 'cpgIslandExt'])
    args = p.parse_args()

    if (args.genome is None or args.bed is None):
        sys.exit(not p.print_help())

    tables = args.tables or ('refGene', 'cpgIslandExt')
    annotate(args.bed, args.genome, tables, in_memory=args.in_memory)

if __name__ == "__main__":
    main()
