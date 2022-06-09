# ------------------------------------------------------------------------------
# Compile ReplicationSuccess_Tutorial
# Leonhard Held
# 2020-01-08
# ------------------------------------------------------------------------------

all: ReplicationSuccess_Tutorial_Brasil.pdf 

ReplicationSuccess_Tutorial_Brasil.pdf: ReplicationSuccess_Tutorial_Brasil.Rnw Makefile
	echo "library('knitr'); knit2pdf('ReplicationSuccess_Tutorial_Brasil.Rnw')" | R --no-save --no-restore
	pdfjam --vanilla --quiet --nup 2x2 --landscape  --outfile ReplicationSuccess_Tutorial_Brasil-nup.pdf ReplicationSuccess_Tutorial_Brasil.pdf



clean:
	rm -f *.dvi *.aux *.log *.toc *.bmt *.bbl *.blg *.idx *.ilg *.ind *.mtc* *.flc *.out *.nav *.vrb *.snm *.gz 
