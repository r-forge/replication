# ------------------------------------------------------------------------------
# Compile ReplicationSuccess_Tutorial
# Leonhard Held
# 2020-01-08
# ------------------------------------------------------------------------------

all: ReplicationSuccess_Tutorial.pdf WorkshopIntro.pdf

ReplicationSuccess_Tutorial.pdf: ReplicationSuccess_Tutorial.Rnw Makefile
	echo "library('knitr'); knit2pdf('ReplicationSuccess_Tutorial.Rnw')" | R --no-save --no-restore
	pdfjam --vanilla --quiet --nup 2x2 --landscape  --outfile ReplicationSuccess_Tutorial-nup.pdf ReplicationSuccess_Tutorial.pdf

WorkshopIntro.pdf: WorkshopIntro.Rnw Makefile
	echo "library('knitr'); knit2pdf('WorkshopIntro.Rnw')" | R --no-save --no-restore
	pdfjam --vanilla --quiet --nup 2x2 --landscape  --outfile WorkshopIntro-nup.pdf WorkshopIntro.pdf


clean:
	rm -f *.dvi *.aux *.log *.toc *.bmt *.bbl *.blg *.idx *.ilg *.ind *.mtc* *.flc *.out *.nav *.vrb
