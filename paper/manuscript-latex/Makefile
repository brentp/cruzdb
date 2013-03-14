cruzdb.pdf: *.tex *.bbl
	rm -f cruzdb.aux cruzdb.pdf cruzdb.log
	pdflatex cruzdb >pdflatex.err
	pdflatex cruzdb >>pdflatex.err

all:
	cruzdb.pdf
