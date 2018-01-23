<!-- README.md is generated from README.Rmd. Please edit that file -->
HierarchicalGOF
---------------

This package presents code and analysis used for assessment of goodness-of-fit for Bayesian hierarchical models. The package is meant to accompany the paper "A Guide to Bayesian Model Checking for Ecologists" by Paul B. Conn, Devin S. Johnson, Perry J. Williams, Sharon R. Melin, and Mevin B. Hooten. Preprint is available at [PeerJ](https://peerj.com/preprints/3390.pdf).

Spatial regression simulations
------------------------------

To run spatial regression simulations, one simply needs to run the script ./inst/SpatialRegression/run\_Poisson\_JAGS\_sims.R. This will iteratively generate data and call JAGS to run code (note you must install JAGS first; see <http://mcmc-jags.sourceforge.net>)

Sea otter N-mixture example
---------------------------

California sea lion example
---------------------------

To run the California sea lion example, the user may use the umbrella function \`run\_attendance\_analysis.' Note this function requires use of RTools to install software from a related github repository (So RTools must be installed).

### Disclaimer

*This software package is developed by scientists at the NOAA Fisheries Alaska Fisheries Science Center and Colorado State University. It should be considered a fundamental research communication. The reccomendations and conclusions presented here are those of the authors and this software should not be construed as official communication by NMFS, NOAA, or the U.S. Dept. of Commerce. In addition, reference to trade names does not imply endorsement by the National Marine Fisheries Service, NOAA. While the best efforts have been made to insure the highest quality, tools such as this are under constant development and are subject to change.*
