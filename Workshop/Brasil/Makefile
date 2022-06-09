# ------------------------------------------------------------------------------
# Compile ReplicationSuccess_Tutorial
# Leonhard Held
# 2020-01-08
# ------------------------------------------------------------------------------

all: Lecture_Brasil.pdf Workshop_Brasil.pdf

Lecture_Brasil.pdf: Lecture_Brasil.Rnw Makefile
	echo "library('knitr'); knit2pdf('Lecture_Brasil.Rnw')" | R --no-save --no-restore
	pdfjam --vanilla --quiet --nup 2x2 --landscape  --outfile Lecture_Brasil-nup.pdf Lecture_Brasil.pdf

Workshop_Brasil.pdf: Workshop_Brasil.Rnw Makefile
	echo "library('knitr'); knit2pdf('Workshop_Brasil.Rnw')" | R --no-save --no-restore
	pdfjam --vanilla --quiet --nup 2x2 --landscape  --outfile Workshop_Brasil-nup.pdf Workshop_Brasil.pdf



clean:
	rm -f *.dvi *.aux *.log *.toc *.bmt *.bbl *.blg *.idx *.ilg *.ind *.mtc* *.flc *.out *.nav *.vrb *.snm *.gz *-concordance.tex
