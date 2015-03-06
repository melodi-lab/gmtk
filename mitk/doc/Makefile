all: mi-toolkit-description.ps mi-toolkit-description.pdf

prev: mi-toolkit-description.dvi
	xdvi mi-toolkit-description.dvi

ps: mi-toolkit-description.ps
pdf: mi-toolkit-description.pdf

mi-toolkit-description.ps: mi-toolkit-description.dvi mi-toolkit-description.bbl
	dvips mi-toolkit-description.dvi -o mi-toolkit-description.ps

mi-toolkit-description.ps.gz: mi-toolkit-description.dvi mi-toolkit-description.bbl
	dvips mi-toolkit-description.dvi -o mi-toolkit-description.ps
	gzip mi-toolkit-description.ps

mi-toolkit-description.pdf: mi-toolkit-description.dvi mi-toolkit-description.bbl
	dvipdf mi-toolkit-description.dvi

mi-toolkit-description.dvi: mi-toolkit-description.tex
	latex mi-toolkit-description

mi-toolkit-description.bbl:
	latex mi-toolkit-description
	bibtex mi-toolkit-description
	latex mi-toolkit-description
	latex mi-toolkit-description

bib:
	latex mi-toolkit-description
	bibtex mi-toolkit-description
	latex mi-toolkit-description
	latex mi-toolkit-description

clean:
	rm -f *~ mi-toolkit-description.log mi-toolkit-description.aux mi-toolkit-description.blg mi-toolkit-description.dvi mi-toolkit-description.ps mi-toolkit-description.ps.gz mi-toolkit-description.pdf

clean-all:
	rm -f *~ mi-toolkit-description.log mi-toolkit-description.aux mi-toolkit-description.blg mi-toolkit-description.dvi mi-toolkit-description.ps mi-toolkit-description.ps.gz mi-toolkit-description.pdf mi-toolkit-description.bbl
