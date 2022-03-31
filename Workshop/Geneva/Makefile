# ------------------------------------------------------------------------------
# Compile ReplicationSuccess_Tutorial
# Leonhard Held
# 2020-01-08
# ------------------------------------------------------------------------------

all: ReplicationSuccess_Geneva.pdf 

ReplicationSuccess_Geneva.pdf: ReplicationSuccess_Geneva.Rnw Makefile
	echo "library('knitr'); knit2pdf('ReplicationSuccess_Geneva.Rnw')" | R --no-save --no-restore
	pdfjam --vanilla --quiet --nup 2x2 --landscape  --outfile ReplicationSuccess_Geneva-nup.pdf ReplicationSuccess_Geneva.pdf



clean:
	rm -f *.dvi *.aux *.log *.toc *.bmt *.bbl *.blg *.idx *.ilg *.ind *.mtc* *.flc *.out *.nav *.vrb *.snm *.gz 
