# Rinbix R implementations and interfaces to In silico research group's 
bioinformatics toolbox [inbix](https://github.com/hexhead/inbix)

Rinbix refers to an R package that ports most of the inbix C++ functions to pure R language functions. The package also includes R interfaces to call the inbix program. The interface accepts and returns R data structures, while taking advantage of C++ performance. (Requires prior installation of the inbix program.) Many support functions from standard bioinformatics packages are included for convenience and to provide a standard calling interface for bioinformatics pipelines/workflows: from data set import and filtering to analysis algorithms to network module interpretation. An example RNA-Seq workflow is included in Rinbix, however, custom workflows can be composed of these convenience functions and “glue” code. As such, Rinbix serves as a testbed for new algorithm and workflow development, and depending on performance needs, serves as a prototype/test implementation for inbix C++ functions. The testbed functions include simulations for synthetic data sets and networks.
 
In most cases, if the resources required are practical, it is best to start exploration of a data set with R, as opposed to starting with C++ right away. However, for large data sets, C++ is often necessary for the filtering/selection stage of the analysis. In this case the R-C++ interface can be used to call the more efficient algorithms while staying in an R workflow start-to-finish.

## Rinbix Build and Install

### CRAN
The CRAN package is under final development. 6/24/17. install.packages(“Rinbix”).

### Github
The latest version of Rinbix is maintained on github. The easiest way to install this version is with devtools::install_github(“insilico/Rinbix”). 

### Source
Download the tarball from github. Install with RStudio (Tools-Install Packages… Install From: Package Archive File) or manually from the command line with: 
 
$ R CMD install Rinbix.xxx.tar.gz
