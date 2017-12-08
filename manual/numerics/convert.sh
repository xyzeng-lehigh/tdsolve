rm numerics.toc
rm numerics.lot
rm numerics.lof
rm numerics.aux

case $1 in
-t) pdflatex numerics
    pdflatex numerics
    ;;
-f) pdflatex numerics
    pdflatex numerics
    bibtex   numerics
    bibtex   numerics
    pdflatex numerics
    pdflatex numerics
    ;;
-p) pdflatex numerics_pdgf
    pdflatex numerics_pdgf
    mv numerics_pdgf.pdf ../
    ;;
*)  pdflatex numerics
    pdflatex numerics
    ;;
esac

mv numerics.pdf ../numerics_opn.pdf
