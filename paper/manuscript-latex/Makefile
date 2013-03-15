cruzdb.pdf: *.tex *.bbl *.bib
	rm -f cruzdb.aux cruzdb.pdf cruzdb.log
	pdflatex cruzdb >pdflatex.err
	bibtex cruzdb.aux
	pdflatex cruzdb >>pdflatex.err

all:
	cruzdb.pdf
