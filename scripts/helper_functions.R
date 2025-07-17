#### Helper functions used in scXpress ####


##################################
####### Loading Queries ##########
##################################

# function for loading .rda ----------------
loadRDa <- function(fileName){
  # Loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
# ------------------------------------------


# ------------------------------------------------------------------------------
# function to load Query (saved as *.rda or *.rds)
# ------------------------------------------------------------------------------
# Prepares counts matrix and if applicable averages expression by specified cluster

load_seurat_query <- function(q_path, avg_by){
  
  ## Load Seurat
  if(grepl("*.rda", q_path, ignore.case = TRUE)){
    q <- loadRDa(q_path)
  }
  if(grepl("*.rds", q_path, ignore.case = TRUE)){
    q <- readRDS(q_path)
  }

  ## Transform Seurat into a dense counts matrix 
  
  ## Preform (if applicable) pre-analysis averaging by cluster
  if (avg_by != "NA"){
    if(avg_by %in% colnames(q@meta.data)) {
      q_dense <- Seurat::AggregateExpression(q, group.by = avg_by) 
      q_dense <- as.matrix(q_dense$RNA,  useNames = FAlSE)
      
      ## Return averaged data
      return(q_dense)
      
    # Warning: could not find avg_by in meta.data  
    } else {
      print(paste("Error: Failed to average expression.", avg_by, "was not in the @meta.data:", q_path))
    }
  } 
  # Return un-averaged data
  q_dense <- as.matrix(Seurat::GetAssayData(q, layer ='counts'), useNames = FAlSE)
  return(q_dense)
}



# ------------------------------------------------------------------------------
# function to load Query (saved as *.csv, *.tsv)
# ------------------------------------------------------------------------------
# Prepares counts matrix and averages expression by specified meta value (if applicable)

load_matrix_query <- function(q_path, m_path, avg_by){
  library(dplyr)
  
  ## Load query 
  q <- data.matrix(as.matrix(data.table::fread(q_path), rownames = 1),rownames.force = TRUE)
  
  ## Load meta
  m <- as.matrix(data.table::fread(m_path), rownames = 1)
  
  ## Preform (if applicable) pre-analysis averaging by cluster
  if (avg_by != "NA") {
    if (avg_by %in% colnames(m)){
      
      # Replace cell barcodes in q with corresponding meta.data column values
      cellIDs <- colnames(q)
      meta_col <- sapply(cellIDs, function(ID)  as.character(m[,avg_by][[ID]]))
      colnames(q) <-meta_col
      
      # Average by column names
      groups <- sort(unique(meta_col))
      
      # Convert to long format for dplyr 
      q_long <- reshape2::melt(q)
      colnames(q_long)<- c("gene", "celltype", "count")
      
      # Summarize values by cell type
      q_aggregated <- q_long %>%
        dplyr::group_by(gene, celltype) %>%
        dplyr::summarize(value = sum(count, na.rm = TRUE)) 
      
      # Convert back to wide format 
      q_wide <- q_aggregated %>%
        tidyr::pivot_wider(names_from = celltype, values_from = value) %>%
        as.matrix()
      
      # Get gene names
      r_names <- q_wide[,1]
      q_wide <- q_wide[,-1]
      
      # transform into numeric 
      q_wide <- apply(q_wide, 2, function(x) as.numeric(x))
      
      # Set genes as row names
      rownames(q_wide) <- r_names
      
      ## Return averaged data
      return(q_wide)
      
    # Warning: could not find avg_by in meta.data   
    } else {
      print(paste("Error: Failed to average expression.", avg_by, "was not in:", m_path))
    }
  }
  
  ## Return un-averaged data
  return(q)
}

# ------------------------------------------------------------------------------
# function to transform build Seurat from counts + meta 
# ------------------------------------------------------------------------------
build_seurat <- function(q){
  ## Transform into Seurat object and normalize raw counts
  q_obj <- Seurat::CreateSeuratObject(counts = q)
  q_obj <- Seurat::NormalizeData(q_obj)
  
  return(q_obj)
}


##################################
####### Loading Signatures #######
##################################

# ------------------------------------------------------------------------------
# function to load External Signatures (using msigdbr):
# ------------------------------------------------------------------------------
load_ex_signatures <- function(species, cat, subcat, gs_style, gene_style){
  
  ## Processing options for msigDB references
  ## Replace NA with NULL 
  if(identical(species, "NA")) species= NULL
  if(identical(cat, "NA")) cat= NULL
  if(identical(subcat, "NA")) subcat= NULL
    
  ## Specifying naming conventions (must match style of sample):
  # gs_style: specifies gene set -> Default: 'gs_name'
  if (identical(gs_style, "NA")) gs_style ="gs_name"
    
  # gene_style: specifies gene name -> Default: 'gene_symbol'
  if (identical(gene_style, "NA")) gene_style ="gene_symbol"
  
  ## Loading reference from msigDB 
  df <- msigdbr::msigdbr(species = species, category = cat, subcategory = subcat)
    
  ## Converting to compatible format
  s <- split(x = df[[gene_style]], f = df[[gs_style]])
  return(s)
}



# ------------------------------------------------------------------------------
# function to isolate specified Signatures:
# ------------------------------------------------------------------------------ 

isolate_signatures <- function(query, ref, gs_name){
  
  ## Check for all_sets option 
  if (identical(gs_name,"all_sets")) {
    gs_list <- ref
  } else {
    # Split names of the Gene Sets
    gs_name <- unlist(strsplit(gs_name,","))
    
    ## Build list of gene sets of interest
    gs_list <- vector(mode = "list", length=0)
    for(i in 1:length(gs_name)){
      name <- gs_name[i]
      
      # Verify that the Gene Set is in the reference
      if(!(name %in% names(ref))){
        print(paste("WARNING:", name, "cannot be found"))
        print(paste(name, "will not be scored"))
      } else {
        # Add the Gene Set to to temp for verification 
        temp <- ref[[name]]
        
        # Verify a proportion of genes are present in the sample (1/4 of signature -> arbitrary)
        missing_genes <- setdiff(temp, rownames(query))
        
        if(length(missing_genes)==length(temp)){
          # Will not score if no genes are present in the sample
          print(paste("WARNING: None of the genes in", name, "are present."))
          print(paste(name, "will not be scored"))
        } else {
          # Will still score if at least one gene is present but prints warning 
          if ((length(missing_genes)/length(temp)) > 0.75){
            print(paste("WARNING: less than 1/4 genes in", name, "are present."))
            print(paste("Verify", name, "is correct"))
          }
          gs_list[[name]] <- temp
        }
      }
    }
  }
  return(gs_list)
}





