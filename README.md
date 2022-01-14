# gRoR
[![GitHub license](https://img.shields.io/github/license/narodkebabci/gRoR?style=plastic)](https://github.com/narodkebabci/gRoR/blob/main/LICENSE)

An enrichment workflow to get rid of redundant mutations

## How it works?

![](https://github.com/narodkebabci/gRoR/blob/media/w.png)
    
## Installation

You can install the latest released version of `gRoR` from Github with:

```R
install.packages("devtools")
devtools::install_github("narodkebabci/gRoR")
```

## Usage

`gRoR` has two main objectives:

1. Construction of class/duplet labels (achieved via `class_formation()` and `duplet_formation()`)
2.1. Compilation of mutation subsets via enrichment method (achieved via `class_reduction()` and `duplet_reduction()`)
2.2. Generation of mutation subsets by undersampling from specified r kcal/mol range (achieved via `ranged_class_reduction()` and `ranged_duplet_reduction()`)


Note that `gRoR` has dependencies, please read DESCRIPTION before using it. 
