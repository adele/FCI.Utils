## FCI.Utils R Package

### Installation

You can download the latest tar.gz file with the source code of the FCI.Utils R package, available at <https://github.com/adele/FCI.Utils/releases/latest>, and install it with the following command, where `path_to_file` represents the full path and file name of the tar.gz file:

``` r
install.packages(path_to_file, repos=NULL, type="source", dependencies=TRUE)
```

Or you can install the development version directly from GitHub. Make sure you have the devtools R package installed. If not, install it with `install.packages("devtools", dependencies=TRUE)`.

``` r
devtools::install_github("adele/FCI.Utils", dependencies=TRUE)
```

Note: if you are asked to update packages, then press "a" for all.

All releases are available at <https://github.com/adele/FCI.Utils/releases/>. If you want a specific version of the FCI.Utils R package, for example, v1.0, you can install it directly from the URL:

First, install the required system libraries. On Debian or Ubuntu install librsvg2-dev:

sudo apt-get install -y librsvg2-dev

or, 


/home/adele/.conda/pkgs/librsvg-2.56.3-h98fae49_0/lib/pkgconfig:/opt/ohpc/pub/mpi/libfabric/1.13.0/lib/pkgconfig:/opt/ohpc/pub/mpi/ucx-ohpc/1.11.2/lib/pkgconfig:/opt/ohpc/pub/mpi/openmpi4-gnu9/4.1.1/lib/pkgconfig

conda update -n base conda
conda update -n causality --all
conda install pkg-config
conda install glib
conda install librsvg
conda install r-xml2
conda install r-desolve
conda install glpk

Then, add pkg-config to PATH and librsvg-2.0.pc to PKG_CONFIG_PATH, e.g.: 

export PKG_CONFIG_PATH=/home/adele/.conda/pkgs/librsvg-2.56.3-h98fae49_0/lib/pkgconfig:$PKG_CONFIG_PATH
export PATH=/home/adele/.conda/pkgs/pkg-config-0.29.2-h36c2ea0_1008/bin/:$PATH


``` r
install.packages('stringi')
install.packages(c("doFuture", "DOT", "rsvg"))
install.packages(c("BFF"))

install.packages(c('lmtest', 'pscl', 'brms', 'MXM'), dependencies=TRUE)

install.packages("https://github.com/adele/FCI.Utils/releases/download/v1.0/FCI.Utils_1.0.tar.gz", repos=NULL, method="libcurl", dependencies=TRUE)
```
