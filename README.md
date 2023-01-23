<img align="right" width="250" src="man/figures/rivnet_logo.png">

# rivnet

Seamless extraction of river networks from Digital Elevation Models data. 

## What this package does 

- Analysis of externally provided DEMs or downloaded from open source repositories (thus interfacing with the `elevatr` package). 
- Extraction of river networks is performed via TauDEM's D8 flow direction algorithm (thus interfacing with the `traudem` package). 
- Resulting river networks are compatible with functions from the `OCNet` package, and can be plotted and analyzed accordingly (and are thus also compatible with packages `igraph` and `SSN`). 
