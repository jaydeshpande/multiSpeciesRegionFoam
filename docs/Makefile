## Makefile — compile physics_guide.tex to PDF
## Requirements: texlive-full (or MiKTeX / MacTeX with siunitx, tcolorbox, booktabs)
##
## Usage:
##   make          — build PDF (runs pdflatex twice for cross-references)
##   make clean    — remove auxiliary files
##   make cleanall — remove auxiliary files AND the PDF

PDF  = physics_guide.pdf
SRC  = physics_guide.tex

.PHONY: all clean cleanall

all: $(PDF)

$(PDF): $(SRC)
	pdflatex -interaction=nonstopmode $(SRC)
	pdflatex -interaction=nonstopmode $(SRC)

clean:
	rm -f *.aux *.log *.toc *.out *.fls *.fdb_latexmk *.synctex.gz

cleanall: clean
	rm -f $(PDF)
