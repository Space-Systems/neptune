.PHONY: neptune.pdf all clean

all: technical-report.pdf

neptune.pdf: technical-report.tex
	latexmk -pdf -pdflatex="pdflatex -shell-escape -interactive=nonstopmode" -use-make neptune.tex
	makeglossaries neptune
	latexmk -pdf -pdflatex="pdflatex -shell-escape -interactive=nonstopmode" -use-make neptune.tex
clean:
	latexmk -CA
