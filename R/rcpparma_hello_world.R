

armamax <- function(A,b=NULL){
	if (is.null(b)) {
		.Call( "armamax", A, PACKAGE = "ArmaUtils" )
	} else {
		.Call( "armamax2", A, b, PACKAGE = "ArmaUtils" )
	}
}

armakron <- function(y,matrices){
	.Call( "armakron", y, matrices, PACKAGE = "ArmaUtils" )
}

