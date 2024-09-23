\name{river-class}

\docType{class}

\alias{river}
\alias{river-class}

\alias{show,river-method}
\alias{$,river-method}
\alias{$<-,river-method}
\alias{names,river-method}
\alias{[[,river,character,missing-method}
\alias{[[<-,river,character,missing-method}
\alias{[[<-,river,character,missing,ANY-method}

\title{river class}

\description{
A \code{river} object contains information on river attributes at different aggregation levels. It can represent a real river network 
(obtained via \code{\link{extract_river}}) or an optimal channel network (obtained via \code{\link[OCNet:create_OCN]{OCNet::create_OCN}}).

The content of a \code{river} object can be treated as a list, hence its objects and sublists can be accessed with both the \code{$} and \code{@} operators. 

 For information on the aggregation levels and on the content of a 
\code{river} object, see \code{\link[OCNet:OCNet-package]{OCNet::OCNet-package}}.
}

\examples{
\dontshow{if (interactive() && traudem::can_register_taudem()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
fp <- system.file("extdata/wigger.tif", package="rivnet")
r <- extract_river(outlet=c(637478,237413),
	DEM=fp)

show(r)
names(r)

# extract or replace parts of a river object
r$dimX
r@dimX
dim <- r[["dimX"]]
r$dimX <- 1
r[["dimX"]]
r[["dimX"]] <- dim

\dontshow{\}) # examplesIf}
}