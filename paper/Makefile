NAME = draft

all: ${NAME}.pdf

${NAME}.aux: ${NAME}.tex
	pdflatex ${NAME}

${NAME}.bbl: ${NAME}.aux
	bibtex ${NAME}

${NAME}.pdf: ${NAME}.tex ${NAME}.bbl
	pdflatex ${NAME}
	pdflatex ${NAME}
