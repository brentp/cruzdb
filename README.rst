A rendered version of the docs is available at: http://pythonhosted.org/cruzdb/

A paper describing cruzdb is in Bioinformatics: http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btt534?ijkey=9I8rQeolKOhzFHv&keytype=ref

cruzdb overview
===============

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
    [(1200999L, 1203486L), (1203543L, 1204010L), (1204082L, 1204420L), (1204682L, 1204836L)]

    # the first 4 exons.
    >>> muc5b.exons[:4]
    [(1200870L, 1200999L), (1203486L, 1203543L), (1204010L, 1204082L), (1204420L, 1204682L)]

    # note that some of these are not coding because they are < cdsStart
    >>> muc5b.cdsStart
    1200929L

    # the extent of the 5' utr.
    >>> muc5b.utr5
    (1200870L, 1200929L)

    # we can get the (first 4) actual CDS's with:
    >>> muc5b.cds[:4]
    [(1200929L, 1200999L), (1203486L, 1203543L), (1204010L, 1204082L), (1204420L, 1204682L)]

    # the cds sequence from the UCSC DAS server as a list with one entry per cds
    >>> muc5b.cds_sequence #doctest: +ELLIPSIS
    ['atgggtgccccgagcgcgtgccggacgctggtgttggctctggcggccatgctcgtggtgccgcaggcag', ...]


    >>> transcript = g.knownGene.filter_by(name="uc001aaa.2").first()
    >>> transcript.is_coding
    False

    # convert a genome coordinate to a local coordinate.
    >>> transcript.localize(transcript.txStart)
    0L

    # or localize to the CDNA position.
    >>> print transcript.localize(transcript.cdsStart, cdna=True)
    None

DataFrames
----------
... are so in. We can get one from a table as::

   >>> df = g.dataframe('knownGene', limit=20) 
   >>> df.columns #doctest: +ELLIPSIS
   Index([name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, proteinID, alignID], dtype=object)



All of the above can be repeated using knownGene annotations by changing 'refGene' to 
'knownGene'. And, it can be done easily for a set of genes.

Spatial
-------

k-nearest neighbors, upstream, and downstream searches are available.
Up and downstream searches use the strand of the query feature to determine the direction:

    >>> nearest = g.knearest("refGene", "chr1", 9444, 9555, k=6)
    >>> up_list = g.upstream("refGene", "chr1", 9444, 9555, k=6)
    >>> down_list = g.downstream("refGene", "chr1", 9444, 9555, k=6)



Mirror
------

The above uses the mysql interface from UCSC. It is now possible to mirror
any tables from UCSC to a local sqlite database via:

   # cleanup

   >>> import os
   >>> if os.path.exists("/tmp/u.db"): os.unlink('/tmp/u.db')

   >>> g = Genome('hg18')



   >>> gs = g.mirror(['chromInfo'], 'sqlite:////tmp/u.db')

and then use as:

   >>> gs.chromInfo
   <class 'cruzdb.sqlsoup.chromInfo'>


Code
----

Most of the per-row features are implemented in `cruzdb/models.py` in the
Feature class. If you want to add something to a feature (like the existing
feature.utr5) add it here.

The tables are reflected using `sqlalchemy`_ and mapped in the
\_\_getattr\_\_\ method of the `Genome` class in `cruzdb/__init__.py`

So a call like::

    genome.knownGene

calls the \_\_getattr\_\_ method with the table arg set to 'knownGene'
that table is then reflected and an object with parent classes of `Feature`
and sqlalchemy's declarative_base is returned.


Contributing
------------

YES PLEASE!

To start coding, it is probably polite to grab your own copy of some of the
UCSC tables so as not to overload the UCSC server. 
You can run something like::

   Genome('hg18').mirror(["refGene", "cpgIslandExt", "chromInfo", "knownGene", "kgXref"], "sqlite:////tmp/hg18.db")

Then the connection would be something like::

    g = Genome("sqlite:////tmp/hg18.db")

If you have a feature you like to use/implement, open a ticket on github for
discussion. Below are some ideas.


.. _`Genomes Database`: http://genome.ucsc.edu/cgi-bin/hgTables
.. _`awk-ful`: https://gist.github.com/1173596
.. _`sqlalchemy`: http://sqlalchemy.org/
