
# Installation Instructions

By Adèle Helena Ribeiro –
[adele.helena\@gmail.com](mailto:adele.helena@gmail.com){.email}


## System Dependencies

### Linux

Before installing the required R packages, ensure that the necessary
system dependencies are installed. Run the following commands in the
shell:

``` sh
sudo apt-get install libmagick++-dev
sudo apt-get install cargo
sudo apt-get install librsvg2-dev
sudo apt-get install libmagick++-dev
sudo apt-get install libavfilter-dev
sudo apt-get install libharfbuzz-dev
sudo apt-get install libcurl4-openssl-dev libssl-dev
sudo apt-get install libudunits2-dev
sudo apt-get install cmake libsodium-dev  
sudo apt-get install libxml2-dev libgdal-dev libfontconfig1-dev libcairo2-dev
```

### macOS

Before starting, create a conda environment and install R (4.3.1) with
the following commands:

``` sh
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -n YOUR_ENV_NAME
conda activate YOUR_ENV_NAME
conda install -c conda-forge r-base=4.3.1
```

Activate the environment and install the necessary packages:

``` sh
conda install -c r r-essentials
conda install pkg-config
conda install glib
conda install -c conda-forge librsvg r-xml2 r-desolve glpk udunits2 r-gsl gsl cxx-compiler r-rsvg
```

Ensure the required system libraries are installed:

``` sh
sudo apt-get install -y librsvg2-dev
```

Then, update `PKG_CONFIG_PATH` and `PATH`:

``` sh
export PKG_CONFIG_PATH=/home/username/.conda/pkgs/librsvg-2.56.3-h98fae49_0/lib/pkgconfig:$PKG_CONFIG_PATH
export PATH=/home/username/.conda/pkgs/pkg-config-0.29.2-h36c2ea0_1008/bin/:$PATH
```

## R Package Installation

Open an R session and execute the following commands:

``` r
install.packages("BiocManager")
BiocManager::install(c("RBGL", "graph", "Rgraphviz"))

package_list <- c('matrixcalc', 'lmtest', 'pscl', 'brms', 'dagitty', 'ggm', 'igraph', 'pcalg', 'SEMgraph', 'doFuture', 'DOT', 'jsonlite', 'rsvg')
install.packages(package_list, dependencies=TRUE, repos='http://cran.us.r-project.org')

new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
if(length(new_packages))
  install.packages(new_packages, dependencies=TRUE, repos='http://cran.us.r-project.org')
```

## Installing BFF version 3.0.1

BFF version 3.0.1 is the last tested and working version. Install it as
follows:

``` sh
wget https://cran.r-project.org/src/contrib/Archive/BFF/BFF_3.0.1.tar.gz
```

In R:

``` r
install.packages(c("BSDA", "hypergeo", "gsl"), dependencies=TRUE)
install.packages("./BFF_3.0.1.tar.gz", repos=NULL, type="source")
```

## Installing MXM version 1.5.5

MXM version 1.5.5 is the last available version. Install it as follows:

``` sh
wget https://cran.r-project.org/src/contrib/Archive/MXM/MXM_1.5.5.tar.gz
```

In R:

``` r
mxm_packages <- c('lme4', 'doParallel', 'relations', 'Rfast', 'visNetwork', 'energy', 'geepack', 'bigmemory', 'coxme', 'Rfast2', 'Hmisc')
mxm_packages <- mxm_packages[!(mxm_packages %in% installed.packages()[,"Package"])]
if(length(mxm_packages))
  install.packages(mxm_packages, dependencies=TRUE, repos='http://cran.us.r-project.org')

install.packages("./MXM_1.5.5.tar.gz", repos=NULL, type="source")
```

## Installing FCI.Utils R Package

Download the latest tar.gz file with the source code of the FCI.Utils R
package, available at
<https://github.com/adele/FCI.Utils/releases/latest>, and install it
with the following command:

``` r
install.packages(path_to_file, repos=NULL, type="source", dependencies=TRUE)
```

Or install the development version directly from GitHub:

``` r
install.packages("devtools", dependencies=TRUE)
devtools::install_github("adele/FCI.Utils", dependencies=TRUE)
```

For a specific version, install it directly from the URL:

``` r
install.packages("https://github.com/adele/FCI.Utils/releases/download/v1.0/FCI.Utils_1.0.tar.gz", repos=NULL, method="libcurl", dependencies=TRUE)
```
