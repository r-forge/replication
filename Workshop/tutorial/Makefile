# ------------------------------------------------------------------------------
# Compile ReplicationSuccess_Tutorial
# Leonhard Held
# 2020-01-08
# ------------------------------------------------------------------------------

build: ReplicationSuccess_Tutorial.Rnw Makefile
	echo "library('knitr'); knit2pdf('ReplicationSuccess_Tutorial.Rnw')" | R --no-save --no-restore

clean:
	rm -f *.dvi *.aux *.log *.toc *.bmt *.bbl *.blg *.idx *.ilg *.ind *.mtc* *.flc *.out *.nav *.vrb
