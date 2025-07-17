#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## Loading Parameters ----------------------------------------------------------

out_file <- args[1]
query_path <- args[2]
signature_path <- args[3]
meta_path <- args[4]
avg_by <- args[5]
gs_name <- args[6] #format: gs_1,gs_2,gs_3

## Loading custom functions ----------------------------------------------------

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

query <- build_seurat(query)

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

## Isolating signatures from the reference 
gs_list <- isolate_signatures(query, ref, gs_name)

## Transform each gene set into a list of genes 
gs_db <- lapply(gs_list, as.list)


## Analysis --------------------------------------------------------------------

## Run Module Score

## Get rid of '-' to follow AddModuleScore naming requirements
old_gs_names <- names(gs_list)
new_gs_names <- unname(sapply(names(gs_list), function(x) gsub("-", ".", x)))

query <- Seurat::AddModuleScore(object = query, features = gs_list, name = paste0(new_gs_names, ".number"), assay = 'RNA', slot = 'data')

# Remove number added by AddModuleScore
newnames <- stringr::str_remove(colnames(query@meta.data), pattern = "\\.number[0-9]+$")
colnames(query@meta.data) <- newnames

## Extract results
q_results <- t(subset(query@meta.data, select = intersect(new_gs_names, (colnames(query@meta.data)))))

## Outputting results ----------------------------------------------------------

## Separating output path and file name 
dir_path <- dirname(out_file)
file_name <- basename(out_file)

## Build path 
dir.create(dir_path, recursive = TRUE)

## Write results to csv
write.csv(q_results, file = out_file, na = "NA", quote = FALSE, row.names = TRUE)
