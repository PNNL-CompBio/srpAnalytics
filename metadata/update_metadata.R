library(tidyverse)
library(xlsx)
library(data.table)
library(webchem)
library(classyfireR)
library(httr)
library(jsonlite)
library(PubChemR)
library(ctxR)
library(rvest)
library(chromote)
library(stringr)

#############################
## PREPARE METADATA SEARCH ##
#############################

# Read all chemical metadata
data = read.xlsx("ChemicalMetadata_Ver_1_1.xlsx", 1) %>%
  mutate(cas = gsub("`", "", cas))

# Identify where we have any chemical data
all_chem_ids = c(
  fread("~/Downloads/all_srp_data/srpCompendiumV5/samplesToChemicals.txt")$Chemical_ID %>% unique(),
  lapply(list.files("~/Downloads/all_srp_data/cleaned_data_v2/", full.names = T), function(x) {
    message(x)
    if (!grepl("Metadata", x)) {return(unique(fread(x)$chemical_id))}
  }) %>% unlist()
) %>% unique()

# Mark where these are in the metadata sheet, as well as where overall "web" column information is missing
data = data %>%
  mutate(
    HaveData = ifelse(chemical_id %in% all_chem_ids, "X", ""), 
    WebClassNA = ifelse(is.na(chem_class_WEB), "X", "")
  ) %>%
  relocate(HaveData, WebClassNA) %>%
  arrange(desc(HaveData), desc(WebClassNA))

# So there is 3 groups of data to consider here:
# 1. Chemicals where we have data, web, and the information just needs to be checked
# 2. Chemicals where we have data, no web, and the information needs to be checked and labels added
# 3. Chemicals where we have no data - we will put these to the side for now
data_web = data %>% filter(HaveData == "X" & WebClassNA == "")
data_noweb = data %>% filter(HaveData == "X" & WebClassNA == "X")
nodata = data %>% filter(HaveData == "" & WebClassNA == "")

############################################
## EXTRACT REQUIRED INFORMATION - COMPTOX ##
############################################

# Pull CAS Numbers where we need more information
cas_list = c(data_web$cas, data_noweb$cas)

# Read through and gather comptox IDs
all_comptox = do.call(bind_rows, lapply(cas_list, function(cas) {
  message(cas)
  tryCatch(ctxR::chemical_starts_with(cas), error = function(e) {NULL})
}))

# Gather information per chemical: dtxsid, dtxcid, casrn, preferredName, smiles
comptox_top_level = all_comptox %>% 
  filter(casrn %in% meta_late$cas) %>% 
  select(dtxsid, dtxcid, casrn, preferredName, smiles) %>% 
  unique() # Only 1289 of the 1345 have matches 

# Gather additional comptox information: inchikey, molFormula, averageMass
comptox_add_info = lapply(gather_by_id$dtxsid, function(x) {
  message(x)
  get_chemical_details(x, Projection = "ccdchemicaldetails") %>%
    dplyr::select(dtxsid, inchikey, molFormula, averageMass)
})

# Finalize comptox information
comptox = do.call(rbind, lapply(comptox_add_info, function(x) {
  x %>% as.character()
})) %>%
  data.frame() %>%
  setNames(c("dtxsid", "inchikey", "molFormula", "averageMass")) %>%
  mutate(averageMass = as.numeric(averageMass)) %>%
  left_join(comptox_top_level)

#########################################################
## EXTRACT REQUIRED INFORMATION - CLASSYFIRE & PUBCHEM ##
#########################################################

# CID Phase---------------------------------------------------------------------

# Pull CIDs
CIDs = lapply(comptox$dtxsid, function(x) {
  message(x)
  vals = PubChemR::get_cids(x)
  if (length(vals$result) > 1) {browser()}
  return(c(x, vals$result[[1]]$result$IdentifierList$CID))
})

# Access cids
cids = do.call(rbind, CIDs) %>%
  data.frame() %>%
  setNames(c("dtxsid", "cid")) %>%
  filter(dtxsid != cid)

# Pull CIDs for convenience
comptox = left_join(comptox, cids, by = "dtxsid")

# Updated tags are really only needed for cases where we do no have labels
comptox_noweb = comptox %>% filter(casrn %in% data_noweb$cas)

# Classyfire phase--------------------------------------------------------------

# Extract all classyfire information
pull_classyfire = function(cid, inchikey) {
  
  # Pull additional identifiers from classyfire
  classyfire = tryCatch({
    
    # Pull the classification information
    cf = classyfireR::get_classification(inchikey)
    cf_df = classyfireR::classification(cf)
    
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
  
  return(c("cid" = as.character(cid), classyfire) %>% enframe() %>% unnest(cols = value) %>% pivot_wider())
  
}

# Run classyfire
classyfire = do.call(bind_rows,
  lapply(1:nrow(rolling), function(el) {
    pull_classyfire(rolling$cid[el], rolling$inchikey[el])
  })
)

# TO DELETE
classyfire = do.call(bind_rows, lapply(list.files("~/Downloads/cid_metadata/", full.names = T), function(x) {
  fread(x) %>% dplyr::select(any_of(c("cid", "superclass", "class", "subclass")))
})) %>%
  filter(!is.na(cid)) %>%
  filter(cid %in% comptox_noweb$cid)

# Roll over
comptox_noweb = left_join(comptox_noweb, classyfire, by = "cid")

# Pubchem tags phase------------------------------------------------------------

#  Extract tags
get_tags <- function(cid) {
  
  # Construct url
  url = paste0("https://pubchem.ncbi.nlm.nih.gov/compound/", cid, "#section=Chemical-Classes")
  
  # Launch chrome and give it a second to read the webpage
  b = ChromoteSession$new()
  b$Page$navigate(url)
  Sys.sleep(5)
  
  # Get the rendered HTML
  html_raw = b$Runtime$evaluate("document.documentElement.outerHTML")$result$value
  b$close()
  
  # Parse with rvest
  page = read_html(html_raw)
  
  # Get all text lines from the page
  all_text = page %>% html_text2()
  lines = unlist(strsplit(all_text, "\n"))
  
  # Filter lines that start with "3.3."
  target_lines = lines[grepl("^\\s*3\\.3\\.", lines)]
  target_lines = trimws(target_lines)
  
  # Extract just the tags (strip the "3.3.X" numbering prefix)
  tags = str_replace(target_lines, "^3\\.3\\.\\d+\\s*", "")
  tags = trimws(tags)
  
  # If tags is longer than 9, stop it
  if (length(tags) > 9) {
    tags = tags[1:9]
  }
  
  message(paste0("cid: ", cid, ", tags: ", paste0(tags, collapse = ", ")))
  return(c("cid" = cid, "tags" = tags))
}

pchem_tags = list(
  list(cid = 329983, tag1 = "Food Additives", tag2 = "Fragrances"),
  list(cid = 62753, tag1 = "Endocrine Disruptors", tag2 = "Flame Retardants", tag3 = "Polymers"),
  list(cid = 6917655, tag1 = "Drugs", tag2 = "Lipids"),
  list(cid = 5315263, tag1 = "Lipids"),
  list(cid = 62074, tag1 = "Drugs", tag2 = "Animal Drugs", tag3 = "Cosmetics", tag4 = "Flavoring Agents", tag5 = "Food Additives", tag6 = "Fragrances"),
  list(cid = 24864132, tag1 = "Pesticides"),
  list(cid = 5280378, tag1 = "Cosmetics", tag2 = "Endocrine Disruptors", tag3 = "Lipids"),
  list(cid = 68827, tag1 = "Drugs", tag2 = "Lipids"),
  list(cid = 107971, tag1 = "Drugs", tag2 = "Endocrine Disruptors", tag3 = "Lipids"),
  list(cid = 6917864, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 392622, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 123596, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 443939, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 16076883, tag1 = "Drugs"),
  list(cid = 6451164, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 66064, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 213039, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 55283, tag1 = "Drugs", tag2 = "Human Drugs", tag3 = "Animal Drugs", tag4 = "Pesticides"),
  list(cid = 64139, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 4463, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 64971, tag1 = "Drugs", tag2 = "Lipids"),
  list(cid = 41684, tag1 = "Drugs", tag2 = "Human Drugs", tag3 = "Animal Drugs"),
  list(cid = 60855, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 5277135, tag1 = "Drugs", tag2 = "Human Drugs", tag3 = "Endocrine Disruptors"),
  list(cid = 3002977, tag1 = "Drugs"),
  list(cid = 60825, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 446541, tag1 = "Drugs", tag2 = "Human Drugs", tag3 = "Endocrine Disruptors"),
  list(cid = 148192, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 92727, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 64142, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 11683, tag1 = "Drugs", tag2 = "Human Drugs", tag3 = "Animal Drugs", tag4 = "Endocrine Disruptors", tag5 = "Lipids"),
  list(cid = 464205, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 479503, tag1 = "Endocrine Disruptors"),
  list(cid = 5281600, tag1 = "Lipids"),
  list(cid = 57469, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 65016, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 37542, tag1 = "Drugs", tag2 = "Human Drugs", tag3 = "Endocrine Disruptors"),
  list(cid = 135398513, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 60877, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 2972, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 60613, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 64927, tag1 = "Drugs", tag2 = "Human Drugs", tag3 = "Animal Drugs"),
  list(cid = 6398764, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 135398741, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 9574768, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 54726191, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 3324, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 58031952, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 60871, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 441386, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 73111, tag1 = "Drugs", tag2 = "Human Drugs", tag3 = "Dietary Ingredients", tag4 = "Lipids"),
  list(cid = 11250133, tag1 = "Lipids"),
  list(cid = 197810, tag1 = "Endocrine Disruptors"),
  list(cid = 135398748, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 5280804, tag1 = "Drugs", tag2 = "Lipids"),
  list(cid = 6451798, tag1 = "Flavoring Agents", tag2 = "Food Additives", tag3 = "Lipids"),
  list(cid = 73659, tag1 = "Lipids"),
  list(cid = 159055, tag1 = "Drugs", tag2 = "Cosmetics", tag3 = "Flavoring Agents", tag4 = "Food Additives", tag5 = "Fragrances"),
  list(cid = 12620, tag1 = "Drugs", tag2 = "Human Drugs", tag3 = "Cosmetics", tag4 = "Endocrine Disruptors", tag5 = "Fragrances", tag6 = "Lipids", tag7 = "Polymers"),
  list(cid = 66828839, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 5311, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 64945, tag1 = "Lipids"),
  list(cid = 164619, tag1 = "Drugs"),
  list(cid = 72281, tag1 = "Drugs", tag2 = "Endocrine Disruptors", tag3 = "Flavoring Agents", tag4 = "Food Additives", tag5 = "Lipids"),
  list(cid = 62881, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 1880, tag1 = "Endocrine Disruptors", tag2 = "Lipids"),
  list(cid = 492405, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 119034, tag1 = "Drugs", tag2 = "Cosmetics"),
  list(cid = 5281607, tag1 = "Endocrine Disruptors", tag2 = "Lipids"),
  list(cid = 25151504, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 10621, tag1 = "Drugs", tag2 = "Cosmetics", tag3 = "Food Additives", tag4 = "Lipids"),
  list(cid = 15165, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 5281605, tag1 = "Endocrine Disruptors", tag2 = "Lipids"),
  list(cid = 4912, tag1 = "Drugs", tag2 = "Human Drugs", tag3 = "Endocrine Disruptors"),
  list(cid = 91656, tag1 = "Endocrine Disruptors", tag2 = "Pesticides"),
  list(cid = 445154, tag1 = "Drugs", tag2 = "Human Drugs", tag3 = "Dietary Ingredients", tag4 = "Endocrine Disruptors", tag5 = "Lipids"),
  list(cid = 45375808, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 1794427, tag1 = "Drugs"),
  list(cid = 5280445, tag1 = "Endocrine Disruptors", tag2 = "Lipids"),
  list(cid = 5092, tag1 = "Endocrine Disruptors"),
  list(cid = 5462355, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 65064, tag1 = "Drugs", tag2 = "Lipids"),
  list(cid = 5281616, tag1 = "Endocrine Disruptors", tag2 = "Lipids"),
  list(cid = 969516, tag1 = "Drugs", tag2 = "Human Drugs", tag3 = "Cosmetics"),
  list(cid = 864, tag1 = "Drugs"),
  list(cid = 2406, tag1 = "Drugs"),
  list(cid = 25154714, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 27944, tag1 = "Drugs", tag2 = "Human Drugs", tag3 = "Animal Drugs"),
  list(cid = 6437380, tag1 = "Drugs", tag2 = "Human Drugs", tag3 = "Endocrine Disruptors"),
  list(cid = 68385, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 90311989, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 9212, tag1 = "Drugs", tag2 = "Endocrine Disruptors"),
  list(cid = 108185, tag1 = "Drugs"),
  list(cid = 65329, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 3034034, tag1 = "Drugs", tag2 = "Human Drugs", tag3 = "Animal Drugs", tag4 = "Cosmetics"),
  list(cid = 736186, tag1 = "Drugs"),
  list(cid = 12947, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 6135, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 68911, tag1 = "Drugs", tag2 = "Human Drugs", tag3 = "Endocrine Disruptors"),
  list(cid = 67505836, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 16129778, tag1 = "Drugs", tag2 = "Human Drugs", tag3 = "Cosmetics", tag4 = "Endocrine Disruptors", tag5 = "Flavoring Agents", tag6 = "Food Additives"),
  list(cid = 6321424, tag1 = "Pesticides"),
  list(cid = 2165, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 135565674, tag1 = "Drugs"),
  list(cid = 216210, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 11556711, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 3000518, tag1 = "Drugs"),
  list(cid = 387447, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 17134, tag1 = "Drugs", tag2 = "Human Drugs", tag3 = "Animal Drugs"),
  list(cid = 64150, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 3957, tag1 = "Drugs", tag2 = "Human Drugs"),
  list(cid = 13258, tag1 = "Polymers"),
  list(cid = 15604, tag1 = "Polymers"),
  list(cid = 87800, tag1 = "Endocrine Disruptors"),
  list(cid = 11796, tag1 = "Endocrine Disruptors")
)


# Acquire pubchem tags
pchem_tags = lapply(comptox_noweb %>% filter(!is.na(cid)) %>% pull(cid), get_tags)

# Add to comptox_noweb
comptox_noweb = left_join(
  comptox_noweb,
  do.call(bind_rows, pchem_tags) %>% setNames(c(
    "cid", "pubchem_class_a", "pubchem_class_b", "pubchem_class_c", "pubchem_class_d",
    "pubchem_class_e", "pubchem_class_f", "pubchem_class_g", 
  ))
)

#######################
## CONSTRUCT OUTPUTS ##
#######################

# Safely updated data.frame without replacing anything that is already defined
safe_update <- function(df, row, colname, value) {
  tryCatch({
    if (is.na(df[row, colname])) {
      df[row, colname] = value
    }
  }, error = function(e) {browser()})
  return(df)
}

# No WEB------------------------------------------------------------------------

# From the data without a web identifier, do the following
# REPLACE: molFormula, averageMass (these are wrong in the original file)
# ADD: cid
# IF NA, then REPLACE: dtxcid, dtxsid, preferredName, smiles, inchikey
#                      superclass, class, subclass, pubmed_class_a-i
data_noweb_final = data_noweb %>%
  select(-c(molFormula, averageMass)) %>%
  left_join(comptox_noweb %>% rename(cas = casrn) %>% select(cas, cid, molFormula, averageMass), by = "cas") %>%
  relocate(cid, .after = cas) %>%
  relocate(averageMass, .after = image_link) %>%
  relocate(molFormula, .after = inchikey)

# Now, let's go line by line and fill in any blanks
for (row in 1:nrow(data_noweb_final)) {
  
  # Grab cas and the metadata
  my_cas = data_noweb_final[row, "cas"]
  my_meta = comptox_noweb[comptox_noweb$cas == my_cas,]
  
  if (nrow(my_meta) == 0) {next}
  
  cols_to_update <- c("dtxcid", "dtxsid", "preferredName", "smiles", "inchikey",
                      "superclass", "class", "subclass",
                      "pubchem_class_a", "pubchem_class_b", "pubchem_class_c", "pubchem_class_d",
                      "pubchem_class_e", "pubchem_class_f", "pubchem_class_g")
  
  for (colname in cols_to_update) {
    data_noweb_final = safe_update(data_noweb_final, row, colname, my_meta[[colname]])
  }
  
}

# WEB---------------------------------------------------------------------------

# From the data without a web identifier, do the following
# REPLACE: molFormula, averageMass (these are wrong in the original file)
# IF NA, then REPLACE: dtxcid, dtxsid, preferredName, smiles, inchikey
data_web_final = data_web %>%
  select(-c(molFormula, averageMass)) %>%
  left_join(comptox %>% rename(cas = casrn) %>% select(cas, cid, molFormula, averageMass), by = "cas") %>%
  relocate(cid, .after = cas) %>%
  relocate(averageMass, .after = image_link) %>%
  relocate(molFormula, .after = inchikey)

# Now, let's go line by line and fill in any blanks
for (row in 1:nrow(data_web_final)) {
  
  # Grab cas and the metadata
  my_cas = data_web_final[row, "cas"]
  my_meta = comptox[comptox$cas == my_cas,]
  
  if (nrow(my_meta) == 0) {next}
  
  cols_to_update <- c("dtxcid", "dtxsid", "preferredName", "smiles", "inchikey")
  
  for (colname in cols_to_update) {
    data_web_final = safe_update(data_web_final, row, colname, my_meta[[colname]])
  }
  
}

## Combine and Write------------------------------------------------------------

data_update_final = bind_rows(
  data_noweb_final,
  data_web_final
)



