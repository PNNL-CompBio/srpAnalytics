library(tidyverse)
library(xlsx)
library(data.table)
library(webchem)
library(classyfireR)
library(httr)
library(jsonlite)
library(PubChemR)

# Read all chemical metadata
data = read.xlsx("ChemicalMetadata_Ver_1_1.xlsx", 1) %>%
  mutate(cas = gsub("`", "", cas))

# Assemble results - function to prevent erroring out
assemble_results <- function(cas = NULL, cid = NULL, preferredName = NULL, chemDescription = NULL,
                             smiles = NULL, inchikey = NULL, molFormula = NULL,
                             averageMass = NULL, superclass = NULL, direct_parent = NULL,
                             class = NULL, subclass = NULL, pubchem_class = NULL) {
  
  base = list("cas" = cas, "cid" = cid, "preferredName" = preferredName, 
               "chemDescription" = chemDescription, "smiles" = smiles, 
               "inchikey" = inchikey, "molFormula" = molFormula, 
               "averageMass" = averageMass, "superclass" = superclass, 
               "direct_parent" = direct_parent, "class" = class, 
               "subclass" = subclass)
  
  if (!is.null(pubchem_class)) {
    tags = pubchem_class
    names(tags) = paste0("pubchem_", letters[seq_along(tags)])
    base = c(base, as.list(tags))
  }
  
  return(base)
}

# Extract missing information
get_chem_info <- function(cas, verbose = TRUE) {
  
  if (verbose) message("Processing CAS: ", cas)
  
  # 1. Get pubchem cid----------------------------------------------------------
  
  # Pull first cid
  cid = tryCatch(get_cid(cas, from = "xref/RN"), error = function(e) NULL)
  
  # If no CID, return all NULL
  if (is.null(cid)) {return(assemble_results())}
  
  # Extract CID and write value
  cid = cid$cid[1]
  if (verbose) message("  PubChem CID: ", cid)
  
  # 2. Get pubchem properties---------------------------------------------------
  
  # Extract CID properties
  props = tryCatch(pc_prop(cid), error = function(e) NULL)
  
  # 3. Extract description------------------------------------------------------
  
  # Extract a description
  description = tryCatch({
    
    url = paste0(
      "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/",
      cid, "/JSON?heading=Record+Description"
    )
    resp = fromJSON(url, flatten = TRUE)
    sections = resp$Record$Section$Section[[1]]$Information
    
    # Extract first available description text
    desc_texts = sections[[1]]$Value.StringWithMarkup
    desc_texts[[1]]$String[1]
    
  }, error = function(e) NULL)
  
  # 4. Classyfire information---------------------------------------------------
  
  # Pull additional identifiers from classyfire
  classyfire = tryCatch({
    
    # Pull the classification information
    cf = get_classification(props$InChIKey)
    cf_df = classification(cf)
    
    get_level = function(level) {
      tryCatch(cf_df$Classification[cf_df$Level == level], error = function(e) return(NULL))
    }
    
    list(
      superclass = get_level("superclass"),
      direct_parent = get_level("direct_parent"),
      class = get_level("class"),
      subclass = get_level("subclass")
    )
    
  }, error = function(e) NULL)
  
  # 5. PubChem classes----------------------------------------------------------
  
  # Extract pubchem information
  pview = get_pug_view(identifier = cid, annotation = "data", domain = "compound")
  pubchem_class = tryCatch(
    {pview$result$Record$Section[[4]]$Section[[4]]$Section %>% lapply(function(x) {x$TOCHeading}) %>% 
        unlist()},
    error = function(e) NULL)
  
  # 6. Return results-----------------------------------------------------------
  
  return(assemble_results(cas = cas,
                          cid = cid,
                          preferredName = props$IUPACName, 
                          chemDescription = description,
                          smiles = props$SMILES,
                          inchikey = props$InChIKey,
                          molForm = props$MolecularFormula,
                          averageMass = props$MolecularWeight,
                          superclass = classyfire$superclass,
                          direct_parent = classyfire$direct_parent,
                          class = classyfire$class,
                          subclass = classyfire$subclass,
                          pubchem_class = pubchem_class))
  
}


for (entry in 1898:nrow(data)) {
  
  cas = data$cas[entry]
  message(entry, "...", cas)
  
  tryCatch({
    get_chem_info(cas, verbose = FALSE) %>%
      enframe() %>%
      unnest(cols = value) %>% 
      pivot_wider() %>%
      fwrite(file.path("~/Downloads/metadata_hunt", paste0(cas, ".txt")), quote = F, row.names = F, sep = "\t")
  }, error = function(e) {message("...not found!")})

}


