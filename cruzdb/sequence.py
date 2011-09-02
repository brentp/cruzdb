import urllib as U

__all__ = ('sequence', )

def _seq_from_xml(xml):
    start = xml.find(">", xml.find("<DNA")) + 1
    end = xml.rfind("</DNA>")
    return xml[start:end].replace(' ', '').replace('\n', '').strip()

def sequence(db, chrom, start, end):
    """
    return the sequence for a region using the UCSC DAS
    server. not the start is 1-based

    >>> sequence('hg18', 'chr2', 2223, 2230)
    'caacttag'
    """
    url = "http://genome.ucsc.edu/cgi-bin/das/%s" % db
    url += "/dna?segment=%s:%i,%i"
    xml = U.urlopen(url % (chrom, start, end)).read()
    return _seq_from_xml(xml)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
