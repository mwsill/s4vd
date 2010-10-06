# TODO: Add comment
# 
# Author: Martin Sill
###############################################################################

# see flexclust package 
#.onLoad <- function(libname, pkgname) {
#	options("s4vd" =
#					list(have_multicore = !inherits(try(loadNamespace("multicore"),
#											silent=TRUE),
#									"try-error")))
#}


setClass('BCs4vd',
		contains = 'BiclustMethod',
		prototype = prototype(
				biclustFunction = function(X,...){s4vd(X,...)}))

BCs4vd <- function() {
	return(new('BCs4vd'))
}


setClass('BCssvd',
		contains = 'BiclustMethod',
		prototype = prototype(
				biclustFunction = function(X,...){ssvdBC(X,...)}))

BCs4vd <- function() {
	return(new('BCssvd'))
}

