all: pubmed.bib
	pdflatex manuscript.tex
	bibtex manuscript.aux
	pdflatex manuscript.tex
	pdflatex manuscript.tex

pubmed.bib: pubmed.lst
	python getref.py pubmed.lst > pubmed.xml
	xsltproc pubmed2bibtex.xsl pubmed.xml > pubmed.bib
