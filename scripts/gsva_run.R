#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## Loading Parameters ----------------------------------------------------------

out_file <- args[1]
query_path <- args[2]
signature_path <- args[3]
meta_path <- args[4]
avg_by <- args[5]
gs_name <- args[6] #format: gs_1,gs_2,gs_3
method <- args[length(args)]

## Loading custom functions ----------------------------------------------------
print(getwd())

source("scripts/helper_functions.R")

## Loading Sample --------------------------------------------------------------

## Sample is a Seurat object (saved as *.rda or *rds):
if(grepl("*.rda", query_path, ignore.case = TRUE)||grepl("*.rds", query_path, ignore.case = TRUE)){
  
  ## Prints warning: all metadata must be pre-included in a Seurat 
  if(meta_path != "NA") print(paste("WARNING: All meta data must be already included in an inputed Seurat"))
  
  ## Load Seurat
  query <- load_seurat_query(query_path, avg_by) 
}


## Sample is a dense counts data table:
if(grepl("*.csv", query_path, ignore.case = TRUE) || grepl("*.tsv", query_path, ignore.case = TRUE)){
  
  ## Load matrix
  query <- load_matrix_query(query_path, meta_path, avg_by)
}


# Loading Lists of Signatures --------------------------------------------------

## External List (using msigdbr):
if(identical(signature_path,"msigdbr")) {
  # Loading Parameters
  species <- args[7]
  category <- args[8]
  subcategory <- args[9]
  gs_style <- args[10]
  gene_style <- args[11]
  
  ## Loading reference
  ref <- load_ex_signatures(species = species, 
                            cat = category, 
                            subcat = subcategory, 
                            gs_style = gs_style, 
                            gene_style = gene_style)

## Internal List (saved as .rds)
} else {
  
  # Loading Parameters
  specifier <- args[7]
  
  ## Loading reference
  if (identical(specifier, "NA")) {
    ref <- readRDS(signature_path)
  } else { 
    ref <- readRDS(signature_path)[[specifier]]
  }
}

## Preprocessing ---------------------------------------------------------------

# Isolating signatures from the reference 
gs_list <- isolate_signatures(query, ref, gs_name)

## Analysis --------------------------------------------------------------------

## Running analysus

if (method == "gsva") q_results <- GSVA::gsva(query, gs_list, method = "gsva", kcdf = "Poisson")
if(method == "gsva_gaussian") q_results <- GSVA::gsva(query, gs_list, method = "gsva", kcdf = "Gaussian")
if (method == "ssgsea") q_results <- GSVA::gsva(query, gs_list, method = "ssgsea", ssgsea.norm = as.logical(FALSE))
if (method == "ssgsea_norm") q_results <- GSVA::gsva(query, gs_list, method = "ssgsea", ssgsea.norm = as.logical(TRUE))

## Outputting results ----------------------------------------------------------

## Separating output path and file name 
dir_path <- dirname(out_file)
file_name <- basename(out_file)

## Build path 
dir.create(dir_path, recursive = TRUE)

## Write results to csv
write.csv(q_results, file = out_file, na = "NA", quote = FALSE, row.names = TRUE)
