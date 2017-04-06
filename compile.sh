#!/bin/bash

echo "Generating Figures, please wait..."
if [[ "$1" == "--figs" ]]; then
    inkscape -D --file=plotgen/structure.svg --export-pdf=plotgen/structure.pdf
    inkscape -D --file=plotgen/flowchart.svg --export-pdf=plotgen/flowchart.pdf
    inkscape -z --file=plotgen/bigiel_extract.svg --export-png=plotgen/bigiel_extract.png --export-dpi=600

    python3 plotgen/toomre_q.py
    python3 plotgen/Q_analysis.py
    python3 plotgen/bigiel_overplot.py
fi

echo " --##--##--##--##--##--##--##--##--##--##--##--##--##--##--##-- "
echo "Running LaTeX"
pdflatex -output-directory=compiled main.tex
bibtex compiled/main.aux
pdflatex -output-directory=compiled main.tex


evince compiled/main.pdf &
