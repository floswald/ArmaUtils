

install:
	R CMD INSTALL ../ArmaUtils --no-multiarch

test:
	Rscript -e "library('testthat'); require(devtools); test('./');"
