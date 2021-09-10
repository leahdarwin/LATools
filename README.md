# LATools

This package is used to construct and analyze locating arrays. A locating array is a combinatorial tool used to design tests that include interactions. This tool set is able to construct locating arrays from a description of test factors. Once the locating array is constructed along and tested to produce results, an algorithm is used to build models. The model builder will inform the user of important factors and interactions in their system.

## Installation

This package is intended to be used in RStudio. To install, first install `devtools`. 

```R
install.packages("devtools")
```

Next, load `devtools` and use it to install this package.

```R
library(devtools)
install_github("/leahdarwin/LATools")
```
## Documentation

All documentation for the package can be accessed in an R interactive session. For example, to access documentation and examples for the function `buildModels()` use the following command.

```R
? buildModels
```

## Authors

This work would not have been possible without the C++ code contribution from Stephen Seidel.
Github: https://github.com/sseidel16/v4-la-tools
