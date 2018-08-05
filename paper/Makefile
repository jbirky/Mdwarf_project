NAME = draft

all: ${NAME}.pdf

${NAME}.pdf: ${NAME}.tex *.bib 
	pdflatex ${NAME}.tex
	rm *.log *.aux *.dvi *.toc *.bbl *.blg