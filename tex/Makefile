RM= /bin/rm -vf

targets= somm.pdf

all: $(targets)

.SECONDARY: $(subst .pdf,.dvi,$(targets))

%.pdf: %.ps
	ps2pdf -dMaxSubsetPct=100 -dCompatibilityLevel=1.2 -dSubsetFonts=true -dEmbedAllFonts=true $<
	$(RM) $<

%.gz: %
	gzip -f $<

%.ps: %.dvi
	dvips -t letter $< -o

%.dvi: %.tex Makefile ann.bib
	latex $<
	latex $<
	bibtex $(subst .tex,,$<)
	latex $<
	latex $<

#table.ps: fitv_parameters.tex header.tex footer.tex
#	cat header.tex fitv_parameters.tex footer.tex > table.tex
#	latex table.tex
#	latex table.tex
#	latex table.tex
#	dvips -Pcmz -o table.ps table.dvi
#	$(RM) table.tex
#
.PHONY: clean spotless

clean:
	$(RM) *.aux *.bbl *.blg *.dvi *.log

spotless: clean
	$(RM) $(targets)
