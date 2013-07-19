cruzdb.pdf: *.tex *.bbl *.bib
	rm -f cruzdb.aux cruzdb.pdf cruzdb.log
	pdflatex cruzdb >pdflatex.err
	bibtex cruzdb.aux
	pdflatex cruzdb >>pdflatex.err

all:
	#cruzdb.pdf
	document.pdf

document.pdf: *.tex *.bbl *.bib
	rm -f document.aux document.pdf document.log document.bbl
	latex document > pdflatex.err
	bibtex document.aux
	latex document >>pdflatex.err
	latex document >>pdflatex.err
	dvipdfm document

pdflatex: *.tex *.bbl *.bib
	rm -f document.aux document.pdf document.log document.bbl
	pdflatex document > pdflatex.err
	bibtex document.aux
	pdflatex document >>pdflatex.err
	pdflatex document >>pdflatex.err
