pdflatex -output-directory=compiled main.tex

inkscape -D --file=plotgen/structure.svg --export-pdf=plotgen/structure.pdf

evince compiled/main.pdf &
