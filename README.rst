cruzdb
======

The UCSC `Genomes Database`_ is a great resource for annoations, regulation
and variation and all kinds of data for a growing number of taxa.
This library aims to make utilizing that data simple so that we can do
sophisticated analyses without resorting to `awk-ful`_, error-prone
manipulations.
As motivation, here's an example of some of the capabilities::

    >>> from cruzdb import Genome

    >>> g = Genome(db="hg18")

    >>> muc5b = g.refGene.filter_by(name2="MUC5B").first()
    >>> muc5b
    refGene(chr11:MUC5B:1200870-1239982)

    >>> muc5b.strand
    '+'

    # the first 4 introns
    >>> muc5b.introns[:4]
    [(1200999, 1203486), (1203543, 1204010), (1204082, 1204420), (1204682, 1204836)]

    # the first 4 exons.
    >>> muc5b.exons[:4]
    [(1200870, 1200999), (1203486, 1203543), (1204010, 1204082), (1204420, 1204682)]

    # note that some of these are not coding because they are < cdsStart
    >>> muc5b.cdsStart
    1200929

    # the extent of the 5' utr.
    >>> muc5b.utr5
    (1200870, 1200929)

    # we can get the (first 4) actual CDS's with:
    >>> muc5b.cds[:4]
    [(1200929, 1200999), (1203486, 1203543), (1204010, 1204082), (1204420, 1204682)]

    # the cds sequence from the UCSC DAS server as a list with one entry per cds
    >>> muc5b.cds_sequence #doctest: +ELLIPSIS
    ['atgggtgccccgagcgcgtgccggacgctggtgttggctctggcggccatgctcgtggtgccgcaggcag', ...]


    >>> transcript = g.knownGene.filter_by(name="uc001aaa.2").first()
    >>> transcript.is_coding
    False

    # convert a genome coordinate to a local coordinate.
    >>> transcript.localize(transcript.txStart)
    0

    # or localize to the CDNA position.
    >>> print transcript.localize(transcript.cdsStart, cdna=True)
    None


this can be repeated using knownGene annotations by changing 'refGene' to 
'knownGene'. And, it can be done easily for a set of genes.

Code
----

Most of the per-row features are implemented in `cruzdb/models.py` in the
Mixin class. If you want to add something to a feature (like the existing
feature.utr5) add it here.

The tables are reflected using `sqlalchemy`_ and mapped in the
\_\_getattr\_\_\ method of the `Genome` class in `cruzdb/__init__.py`

So a call like::

    genome.knownGene

calls the \_\_getattr\_\_ method with the table arg set to 'knownGene'
that table is then reflected and an object with parent classes of `Mixin`
and sqlalchemy's declarative_base is returned.


Contributing
------------

YES PLEASE!

To start coding, it is probably polite to grab your own copy of some of the
UCSC tables so as not to overload the UCSC server. After setting up mysql,
you can run this `script <https://gist.github.com/987144>`_ to get your
own set of commonly used tables.

Then the connection line above would be something like::

    g = Genome(db="hg18", host="localhost", user="brentp")

If you have a feature you like to use/implement, open a ticket on github for
discussion. Below are some ideas.


TODO
----

 + easily make a local copy of a set of tables--to, e.g. sqlite
   Genome.mirror('hg18', ['refGene', 'knownwGene', 'cpgIslantExt'], to_url='sqlite:///hg18.db')

 + examples / THINGS this should make easy:
 + https://lists.soe.ucsc.edu/pipermail/genome/2011-August/026941.html
 + https://lists.soe.ucsc.edu/pipermail/genome/2011-August/026899.html
 + http://biostar.stackexchange.com/questions/10144/ucsc-mysql-query/10147#10147
 + https://lists.soe.ucsc.edu/pipermail/genome/2011-September/027053.html
 + simple name lookup/conversion ens -> hugo, etc.

 + methods for nearest gene
 + use bin to do more efficent query
 + auto join to kgXref
 + histone, cpg, other informative tracks
 + intersect annos with snps

 + sphinx docs
 + **tests tests tests**
 + useful scripts in scripts/ dir


.. _`Genomes Database`: http://genome.ucsc.edu/cgi-bin/hgTables
.. _`awk-ful`: https://gist.github.com/1173596
.. _`sqlalchemy`: http://sqlalchemy.org/
