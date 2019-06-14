################################################################################
## Build/Check/Install a Package
## Manuela Ott, adapted from Sebastian Meyer
## 2018-11-21
################################################################################
     
R := R
PKG := ReplicationSuccess
VERSION := $(strip $(shell grep "^Version:" pkg/DESCRIPTION | cut -f 2 -d ":"))

build:
	$R CMD build pkg

check: build
	_R_CHECK_FORCE_SUGGESTS_=FALSE $R CMD check ${PKG}_${VERSION}.tar.gz
	@cd ${PKG}.Rcheck; nwarn=`grep -c "^Warning" ${PKG}-Ex.Rout`; \
	if [ $$nwarn -gt 0 ]; then echo "\n\tWARNING: $$nwarn" \
        "warning(s) thrown when running examples,\n" \
	"\t         see file ${PKG}.Rcheck/${PKG}-Ex.Rout\n"; fi

cran: build
	_R_CHECK_FORCE_SUGGESTS_=FALSE $R CMD check --as-cran ${PKG}_${VERSION}.tar.gz

install: build
	$R CMD INSTALL ${PKG}_${VERSION}.tar.gz

#manual:
#	$R CMD Rd2pdf --batch --force --output=${PKG}.pdf ${PKG}

## all targets are "phony"
.PHONY: build check install
#.PHONY: build check install manual
