  setwd("D:/Landscape_210/candidate gene")
  file <- "cicar.CDCFrontier.gnm3.ann1.NPD7.gcv_genes.gff3 (1)-LAPTOP-S8OS68Q5.xlsx"
  library(readr)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(lme4)
  library(tibble)
  library(readxl)
  # Open the main sheet
  genes <- read_excel(
    file,
    sheet = "cicar.CDCFrontier.gnm3.ann1.NPD"
  )
  library(dplyr)
  library(purrr)
  library(stringr)
  
  genes_compact <- genes %>%
    filter(type == "gene") %>%
    mutate(
      # Reconstruct all annotation fields into one string
      annotation_full = pmap_chr(
        select(., attributes, starts_with("...")),
        ~ paste(
          na.omit(c(...)),
          collapse = ";"
        )
      ),
      
      # Extract the gene Name
      gene_name = str_extract(
        annotation_full,
        "(?<=Name=)[^;]+"
      ),
      
      # Extract all GO terms
      GO = str_extract_all(
        annotation_full,
        "GO:\\d{7}"
      ) %>%
        map(~ paste(unique(.x), collapse = ";"))
    ) %>%
    select(
      chromosome,
      chr,
      source,
      type,
      start,
      stop,
      score,
      strand,
      phase,
      gene_name,
      GO,
      annotation_full
    )
  

  
  
  library(dplyr)
  library(stringr)
  
  # Candidate positions
  positions <- c(
    "S6_3841200",
    "S2_47682119",
    "S3_4383165",
    "S6_12411634",
    "S7_17650643"
  )
  
  # Convert positions to chromosome + bp position
  regions <- tibble(
    marker = positions,
    chr = str_extract(positions, "(?<=S)\\d+"),
    position = as.numeric(str_extract(positions, "(?<=_)\\d+"))
  ) %>%
    mutate(
      region_start = pmax(0, position - 500000),
      region_end   = position + 500000
    )
  
  regions
  
  genes_nearby <- genes_compact %>%
    mutate(
      chr_num = str_extract(chr, "\\d+")
    ) %>%
    inner_join(
      regions,
      by = c("chr_num" = "chr"),
      relationship = "many-to-many"
    ) %>%
    filter(
      start <= region_end,
      stop >= region_start
    ) %>%
    select(
      marker,
      position,
      region_start,
      region_end,
      gene_name,
      GO,
      chr,
      start,
      stop,
      strand,
      everything()
    )
  
  
  
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  genes_GO_long <- genes_nearby %>%
    select(
      marker,
      gene_name,
      chr,
      start,
      stop,
      strand,
      GO
    ) %>%
    unnest_longer(GO) %>%
    separate_rows(GO, sep = ";") %>%
    filter(
      !is.na(GO),
      GO != "",
      str_detect(GO, "^GO:\\d{7}$")
    ) %>%
    distinct(
      marker,
      gene_name,
      GO
    )
  
  head(genes_GO_long)  
  
  
  
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  
  background_GO <- genes_compact %>%
    select(
      gene_name,
      GO
    ) %>%
    separate_rows(GO, sep = ";") %>%
    filter(
      !is.na(gene_name),
      gene_name != "",
      !is.na(GO),
      GO != "",
      str_detect(GO, "^GO:\\d{7}$")
    ) %>%
    distinct(
      gene_name,
      GO
    )
  
  # Check the background
  background_GO %>%
    summarise(
      n_genes = n_distinct(gene_name),
      n_GO_terms = n_distinct(GO)
    )

  
  candidate_genes <- genes_GO_long %>%
    distinct(gene_name) %>%
    pull(gene_name)
  
  length(candidate_genes)  

  BiocManager::install(
    "enrichplot",
    dependencies = TRUE,
    ask = FALSE,
    update = FALSE,
    force = TRUE
  )
  library(enrichplot)
  library(clusterProfiler)
 

    
  
  library(GO.db)
  library(AnnotationDbi)
  library(dplyr)
  library(stringr)
  
  # ---------------------------------------------------------
  # 1. Retrieve GO descriptions for all candidate annotations
  # ---------------------------------------------------------
  
  go_info <- AnnotationDbi::select(
    GO.db,
    keys = unique(genes_GO_long$GO),
    keytype = "GOID",
    columns = c("TERM", "ONTOLOGY")
  ) %>%
    distinct(GOID, .keep_all = TRUE)
  
  
  # ---------------------------------------------------------
  # 2. Add GO description and ontology to candidate genes
  # ---------------------------------------------------------
  
  candidate_GO <- genes_GO_long %>%
    left_join(
      go_info,
      by = c("GO" = "GOID")
    )
  
  
  # ---------------------------------------------------------
  # 3. Keywords related to environmental adaptation
  # ---------------------------------------------------------
  
  adaptation_keywords <- paste(
    c(
      # General adaptation / environment
      "adapt",
      "stress",
      "abiotic",
      "biotic",
      "environment",
      "stimulus",
      
      # Temperature
      "temperature",
      "heat",
      "cold",
      
      # Water
      "drought",
      "water deprivation",
      "water stress",
      "osmotic",
      
      # Salinity
      "salt",
      "salinity",
      
      # Light / photoperiod
      "light",
      "photoperiod",
      "circadian",
      "rhyth",
      "season",
      
      # Phenology
      "flower",
      "floral",
      "flowering",
      
      # Hormones
      "hormone",
      "auxin",
      "abscisic",
      "ethylene",
      "jasmon",
      "salicylic",
      
      # Oxidative stress
      "oxidative",
      "reactive oxygen",
      "redox",
      
      # Defence
      "defense",
      "defence",
      
      # Signalling
      "signal",
      "signaling",
      "signalling",
      
      # Resource acquisition / transport
      "transport",
      "nutrient",
      "nitrogen",
      "phosphate",
      "ion"
    ),
    collapse = "|"
  )
  
  
  # ---------------------------------------------------------
  # 4. Identify adaptation-related GO annotations
  #    WITHOUT restricting ontology
  # ---------------------------------------------------------
  
  adaptation_GO <- candidate_GO %>%
    filter(
      !is.na(TERM),
      str_detect(
        str_to_lower(TERM),
        adaptation_keywords
      )
    ) %>%
    distinct(
      marker,
      gene_name,
      GO,
      TERM,
      ONTOLOGY
    ) %>%
    arrange(
      marker,
      gene_name,
      ONTOLOGY,
      GO
    )
  
  adaptation_gene_summary <- adaptation_GO %>%
    group_by(marker, gene_name) %>%
    summarise(
      adaptation_GO = paste(unique(GO), collapse = "; "),
      adaptation_terms = paste(unique(TERM), collapse = "; "),
      ontology = paste(
        unique(ONTOLOGY),
        collapse = "; "
      ),
      .groups = "drop"
    )

  
  adaptation_GO_classified <- candidate_GO %>%
    mutate(
      term_lower = str_to_lower(TERM),
      
      adaptation_category = case_when(
        
        # --------------------------------
        # VERY STRONG / DIRECT ADAPTATION
        # --------------------------------
        
        str_detect(
          term_lower,
          "temperature|heat|cold"
        ) ~ "Temperature response",
        
        str_detect(
          term_lower,
          "drought|water deprivation|water stress|water deficit"
        ) ~ "Drought / water stress",
        
        str_detect(
          term_lower,
          "osmotic"
        ) ~ "Osmotic stress",
        
        str_detect(
          term_lower,
          "salt stress|salinity|response to salt"
        ) ~ "Salinity stress",
        
        str_detect(
          term_lower,
          "photoperiod|photoperiodic"
        ) ~ "Photoperiod response",
        
        str_detect(
          term_lower,
          "response to light|light response"
        ) ~ "Light response",
        
        str_detect(
          term_lower,
          "circadian|circadian rhythm"
        ) ~ "Circadian regulation",
        
        str_detect(
          term_lower,
          "flowering|flower development|floral"
        ) ~ "Flowering / phenology",
        
        str_detect(
          term_lower,
          "response to hormone|hormone response"
        ) ~ "Hormone response",
        
        str_detect(
          term_lower,
          "defense response|defence response"
        ) ~ "Defense response",
        
        # --------------------------------
        # OXIDATIVE / REDOX ADAPTATION
        # --------------------------------
        
        str_detect(
          term_lower,
          "oxidative stress|reactive oxygen|reactive nitrogen"
        ) ~ "Oxidative stress",
        
        str_detect(
          term_lower,
          "redox homeostasis|cell redox homeostasis"
        ) ~ "Redox homeostasis",
        
        str_detect(
          term_lower,
          "peroxiredoxin|thioredoxin"
        ) ~ "Redox protection",
        
        str_detect(
          term_lower,
          "oxidoreductase activity"
        ) ~ "Oxidoreductase activity",
        
        # --------------------------------
        # ION / ELECTROPHYSIOLOGICAL
        # --------------------------------
        
        str_detect(
          term_lower,
          "potassium ion transport|potassium ion transmembrane"
        ) ~ "Potassium transport",
        
        str_detect(
          term_lower,
          "calcium ion transport|calcium signaling|calcium-mediated"
        ) ~ "Calcium signaling / transport",
        
        str_detect(
          term_lower,
          "ion homeostasis"
        ) ~ "Ion homeostasis",
        
        # --------------------------------
        # SIGNALING / REGULATION
        # --------------------------------
        
        str_detect(
          term_lower,
          "small gtpase-mediated signal transduction"
        ) ~ "Small GTPase signaling",
        
        str_detect(
          term_lower,
          "protein phosphorylation|protein dephosphorylation"
        ) ~ "Protein phosphorylation signaling",
        
        str_detect(
          term_lower,
          "ubiquitination|deubiquitination"
        ) ~ "Protein ubiquitination",
        
        # --------------------------------
        # CELLULAR / STRUCTURAL RESPONSE
        # --------------------------------
        
        str_detect(
          term_lower,
          "cell wall modification"
        ) ~ "Cell wall remodeling",
        
        str_detect(
          term_lower,
          "actin cytoskeleton organization"
        ) ~ "Cytoskeleton remodeling",
        
        # --------------------------------
        # OTHER
        # --------------------------------
        
        TRUE ~ NA_character_
      )
    )  

  
  adaptation_high_confidence <- adaptation_GO_classified %>%
    filter(
      adaptation_category %in% c(
        "Temperature response",
        "Drought / water stress",
        "Osmotic stress",
        "Salinity stress",
        "Photoperiod response",
        "Light response",
        "Circadian regulation",
        "Flowering / phenology",
        "Hormone response",
        "Defense response",
        "Oxidative stress",
        "Redox homeostasis",
        "Redox protection",
        "Oxidoreductase activity"
      )
    ) %>%
    distinct(
      marker,
      gene_name,
      GO,
      TERM,
      ONTOLOGY,
      adaptation_category
    ) 
  