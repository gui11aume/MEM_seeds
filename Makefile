all: pubmed.bib
	pdflatex MEM.tex
	bibtex MEM.aux
	pdflatex MEM.tex
	pdflatex MEM.tex

pubmed.bib: pubmed.lst
	python getref.py pubmed.lst > pubmed.xml
	xsltproc pubmed2bibtex.xsl pubmed.xml > pubmed.bib
