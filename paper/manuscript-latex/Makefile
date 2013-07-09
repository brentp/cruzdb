cruzdb.pdf: *.tex *.bbl *.bib
	rm -f cruzdb.aux cruzdb.pdf cruzdb.log
	pdflatex cruzdb >pdflatex.err
	bibtex cruzdb.aux
	pdflatex cruzdb >>pdflatex.err

all:
	#cruzdb.pdf
	document.pdf

document.pdf: *.tex *.bbl *.bib
	rm -f document.aux document.pdf document.log
	pdflatex document > pdflatex.err
	bibtex document.aux
	pdflatex document >>pdflatex.err
