# Usage: sh convert.sh [name]      (name of the tex file without .tex)
#   sh convert.sh             builds manual.tex    -> ../manual.pdf
#   sh convert.sh potassium   builds potassium.tex -> ../potassium.pdf
name=${1:-manual}
pdflatex $name
pdflatex $name
mv $name.pdf ../
