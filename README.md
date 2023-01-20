# rivers

Seamless extraction of river networks from Digital Elevation Models data. 

## What this package does 

- Analysis of externally provided DEMs or downloaded from open source repositories (thus interfacing with the \code{elevatr} package). 
- Extraction of river networks is performed via TauDEM's D8 flow direction algorithm (thus interfacing with the \code{traudem} package). 
- Resulting river networks are compatible with functions from the 'OCNet' package, and can be plotted and analyzed accordingly. 