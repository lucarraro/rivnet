\name{get_river_data}
\alias{get_river_data}
\title{Get river data}
\usage{
get_river_data(outlet=c(735010, 261530),
                           ext=c(700000, 770000, 220000, 270000),
                           EPSG=21781,
                           zoom=9,
						   as.channelNetwork=TRUE,
						   as.rast=FALSE,
						   filename=NULL)
}
\description{
Downloads a river
}
