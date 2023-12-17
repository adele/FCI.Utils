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

All releases are available at <https://github.com/adele/FCI.Utils/releases/>. If you want a specific version of the FCI.Utils R package, for example, v1.1, you can install it directly from the URL:

``` r
install.packages("https://github.com/adele/FCI.Utils/releases/download/v1.0/FCI.Utils_1.1.tar.gz", repos=NULL, method="libcurl", dependencies=TRUE)
```
