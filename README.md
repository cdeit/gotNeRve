# gotNeRve

`gotNeRve` is an R package implementing a gene expression-based classifier for predicting high or low **s**tromal **t**umor-**a**ssociated **n**erve (**STaN**) abundance in non-muscle-invasive bladder cancer (NMIBC).

The classifier uses a 12-gene transcriptomic signature derived from bulk RNA-sequencing data and can be applied to individual samples or larger cohorts.

## Installation

You can install the development version of `gotNeRve` from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("YOURUSERNAME/gotNeRve")
```

## Usage

`gotNeRve()` accepts gene expression matrices with genes as rows and samples as columns. Raw counts, CPM-normalized expression, and log2-CPM expression are supported.

```r
library(gotNeRve)

# Raw counts
pred <- gotNeRve(my_counts, input_type = "counts")

# CPM
pred <- gotNeRve(my_cpm, input_type = "cpm")

# log2-CPM
pred <- gotNeRve(my_log2cpm, input_type = "log2cpm")
```

The function returns the predicted STaN class and class probabilities for each sample.

Raw counts are recommended when available because `gotNeRve()` performs the preprocessing required by the classifier internally.

## Model

The bundled classifier is a locked 12-gene elastic-net model developed in an NMIBC RNA-sequencing cohort and evaluated in an independent RNA-sequencing cohort with matched histologic STaN quantification.

The complete model-development and evaluation workflow is provided in the package vignette.

## Citation

If you use `gotNeRve` in your research, please cite the following publications:

> Deiter CS, de Jong FC, Olislagers M, Jordan KR, Zuiverloon TCM, Costello JC.  
> _Stromal innervation is prognostic in non-muscle-invasive bladder cancer and predicted by a distinct microenvironmental transcriptomic signature._  
> bioRxiv. 2026. [Citation/DOI to be added]

> de Jong FC, Laajala TD, Hoedemaeker RF, Jordan KR, van der Made ACJ, Boevé ER, van der Schoot DKE, Nieuwkamer B, Janssen EAM, Mahmoudi T, Boormans JL, Theodorescu D, Costello JC, Zuiverloon TCM. 
> _Non-muscle-invasive bladder cancer molecular subtypes predict differential response to intravesical Bacillus Calmette-Guérin._ Sci Transl Med. 2023 May 24;15(697):eabn4118. 
> doi: [10.1126/scitranslmed.abn4118](https://doi.org/10.1126/scitranslmed.abn4118)

## License

[Add license information]