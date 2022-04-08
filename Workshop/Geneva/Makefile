# ------------------------------------------------------------------------------
# Compile ReplicationSuccess_Tutorial
# Leonhard Held
# 2020-01-08
# ------------------------------------------------------------------------------

all: Lecture_Geneva.pdf Workshop_Geneva.pdf

Lecture_Geneva.pdf: Lecture_Geneva.Rnw Makefile
	echo "library('knitr'); knit2pdf('Lecture_Geneva.Rnw')" | R --no-save --no-restore
	pdfjam --vanilla --quiet --nup 2x2 --landscape  --outfile Lecture_Geneva-nup.pdf Lecture_Geneva.pdf

Workshop_Geneva.pdf: Workshop_Geneva.Rnw Makefile
	echo "library('knitr'); knit2pdf('Workshop_Geneva.Rnw')" | R --no-save --no-restore
	pdfjam --vanilla --quiet --nup 2x2 --landscape  --outfile Workshop_Geneva-nup.pdf Workshop_Geneva.pdf



clean:
	rm -f *.dvi *.aux *.log *.toc *.bmt *.bbl *.blg *.idx *.ilg *.ind *.mtc* *.flc *.out *.nav *.vrb *.snm *.gz *-concordance.tex
