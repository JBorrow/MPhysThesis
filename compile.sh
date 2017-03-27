pdflatex -output-directory=compiled main.tex

inkscape -D --file=plotgen/structure.svg --export-pdf=plotgen/structure.pdf
inkscape -D --file=plotgen/flowchart.svg --export-pdf=plotgen/flowchart.pdf

evince compiled/main.pdf &
