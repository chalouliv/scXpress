------------------------------------------------------------------------

# scXpress

scXpress is a snakemake pipeline designed to analyze gene set (signature) enrichment in single-cell RNA sequencing data. It supports three methods of enrichment scoring: single-sample Gene Set Enrichment Analysis (ssGSEA), Gene Set Variation Analysis (GSVA), and Module Score.

## Set-up:

1.  [Workflow](#1-Workflow)
2.  [Prepare signatures](#2-Prepare-signatures)
3.  [Prepare queries](#3-Prepare-queries)
4.  [Decide methods](#4-Decide-methods)
5.  [Prepare configfile](#5-Prepare-configfile)
6.  [Installation and dependencies](#6-Installation-and-dependencies)
7.  [Prepare HPC submission script](#7-Prepare-HPC-submission-script)
8.  [Tips](#8-Tips)

## 1. Workflow

scXpress is designed to run each specified method for all queries and for all signatures.

<img src="https://github.com/chalouliv/scXpress/blob/main/graphics/scXpress_workflow.png" width="600"/>

## 2. Prepare signatures

Ultimately signatures must be a list of gene names. There are 2 options available to achieve this:

### Internal List of Signatures

To reference a custom list of signatures, input a path to the list in the configfile (.rds).

```path: /path/to/list.rds```

This list must be formatted as a list of character vectors, where each vector is named after the signature. The vector itself contains the genes associated with that signature.

### MsigDB

To pull signatures from the [MsigDB](https://www.gsea-msigdb.org/gsea/msigdb), input "msigdbr" in place of a path.

```path: msigdbr```

From there specify species, category, sub_category, gene nomenclature, and signature nomenclature as in the [msigdbr()](https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html) function.

## 3. Prepare queries

A query may be imputed as a **Seurat object** or a **counts matrix**.

The naming convention for genes in the query must match the naming convention in the signatures.

**Seurat object** (.rds or .rda)

A Seurat object must have a RNA `counts` layer saved in `assays`.

`Meta.data:` Any meta.data must already be included in the Seurat object.

**Counts matrix** (.csv or .tsv)

A matrix must be a *cell x gene matrix* of RNA counts.

``` bash
 '', gene1, gene2, gene3
 cell1, - , - , -
 cell2, - , - , -
 cell3, - , - , -
 cell4, - , - , -
```

`Meta.data:` There is the option to include a separate meta.data file which can be used for pre- or post-averaging. This file should be a .csv or .tsv with the following format.

``` bash
 '', meta_category1, meta_category2, meta_category3
 cell1, - , - , -
 cell2, - , - , -
 cell3, - , - , -
 cell4, - , - , -
```

### Pre-analysis clustering

There is an option to sum expression across predefined groups (clusters) before performing enrichment scoring. These groups must be defined by a column in the ```meta.data``` and specified using ```avg_by``` in the [configfile](#5-Prepare-configfile).

## 4. Decide methods

scXpress has currently 3 methods available for calculating a gene signature's enrichment score:

### Single-Sample Gene Set Enrichment Analysis (ssGSEA)

<img src="https://github.com/chalouliv/scXpress/blob/main/graphics/ssgsea.png" width="600"/>

**Normalization:** By default, this method does not include the final normalization step outlined in Barbie et al, Nature (2009). This allows scores to be calculated in each cell entirely independently from the rest of the query. There is an option to include this final step.

*From: Barbie et al, Nature (2009)*

### Gene Set Variation Analysis (GSVA)

<img src="https://github.com/chalouliv/scXpress/blob/main/graphics/gsva.png" width="600"/>

**Distribution:**  The first step of GSVA normalizes gene expression across the query using an expression statistic calculated with a Kernel Cumulative Distribution Function (KCDF). Understanding the query data is crucial, as it determines which distribution should be used. By default, Poisson distribution is used, as it is appropriate for integer count data. For continuous expression data, Gaussian distribution should be chosen instead.

*From: Hänzelmann, Castelo, and Guinney, BMC Bioinformatics (2013)*

### Module Score

<img src="https://github.com/chalouliv/scXpress/blob/main/graphics/mod_score.png" width="600"/>

*From: Tirosh et al, Science (2016)*\`

## 5. Prepare configfile

See: [Example Config](examples/config.example)

## 6. Installation and dependencies

There are two options when using scXpress.

### With an Apptainer image (suggested)

scXpress is designed to run within a container created using an Apptainer (v.1.2.4) image, which is available on Docker Hub at [chalouliv/scxpress:test_v6](https://hub.docker.com/layers/chalouliv/scxpress/test_v6/images/sha256-83b7e4cf133a492ff4d22ba34afbfea665fbbd7ae0440b85725bf5843853d3f2?context=repo). This container is a custom environment that has all necessary dependencies.

To run scXpress with the container, first build the image from Docker Hub:

``` bash
apptainer build <image name> docker://chalouliv/scxpress:test_v6
```

```NOTE:``` If you are apart of the Kleinman group, the image has been pre-built at the following locations:

```         
Hydra: /project/kleinman/charlotte.livingston/from_hydra/scxpress.sif

Narval: /lustre06/project/6004736/chaliv/from_narval/scxpress.sif
```

Next see [Prepare HPC submission script](#7-Prepare-HPC-submission-script) to run scXpress.

### With Manual Installation

It is also possible to install the dependencies of scXpress manually.

scXpress has been designed with R version 4.2.2 and Python version 3.11.2.

**Python modules**

``` bash
pip install snakemake # v.8.16.0
```

**R packages**

Installing Seurat:

``` bash
## Seurat requires an older version of Matrix
install.packages("remotes") # v.2.4.2
remotes::install_version("Matrix", version = "1.6.5", repos = "https://cran.r-project.org")

install.packages("Seurat") # v.4.3.0 (SeuratObject v4.1.3)
```

Installing GSVA:

``` bash
install.packages("BiocManager") # v.1.30.20
BiocManager::install("GSVA") # v.1.46.0
```

Other packages:

``` bash
pkgs <- c("data.table", # v.1.14.8
          "stringr", # v.1.5.0
          "dplyr", # v.1.1.0  
          "ggplot2", #v.3.4.1  
          "reshape2", #v.1.4.4 
          "RColorBrewer", #v.1.1-3
          "msigdbr") #v.7.5.1
install.packages(pkgs)
```

Next see [Prepare HPC submission script](#7-Prepare-HPC-submission-script) to run scXpress.

## 7. Prepare HPC submission script

With the [Apptainer image](#6-Installation-and-dependencies), See: [submission_w.sh](examples/submission_w_container.sh)

With [manual installation](#6-Installation-and-dependencie), See: [submission_wo.sh](examples/submission_wo_container.sh)

## 8. Tips

**Apptainer**

To test the pipeline before submitting a job, run an interactive Apptainer shell

``` bash
srun --nodes=1 --mem=20G --pty bash -i # Request an interactive node from HPC
apptainer shell <image_file> 
```

**Snakemake**

Before submitting, execute a dry run

``` bash
snakemake -s <snakefile> --configfile <configfile> -n
```

