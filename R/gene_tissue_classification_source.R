
# Source code for 
#http://www.bioconductor.org/packages/devel/bioc/vignettes/HPAanalyze/inst/doc/b_HPAanalyze_indepth.html

if (!requireNamespace("HPAanalyze", quietly = TRUE)) {
  BiocManager::install("HPAanalyze")
}
library(HPAanalyze)
library(tibble)

downloadedData <- hpaDownload(downloadList='histology')

# # Define the categories and tissues as a tibble, including prostate and spleen
# tissues_tbl <- tibble::tibble(
#   Category = c("Respiratory System", "Digestive System", "Reproductive System", "Urinary System", 
#                "Endocrine System", "Integumentary System", "Musculoskeletal System", 
#                "Nervous System", "Lymphatic System"),
#   Tissues = c("bronchus, lung, nasopharynx, esophagus",
#               "duodenum, rectum, small intestine, colon, pancreas, liver, gallbladder, appendix, stomach",
#               "fallopian tube, testis, ovary, cervix, endometrium, vagina, prostate, epididymis, seminal vesicle, breast, placenta",
#               "urinary bladder, kidney",
#               "adrenal gland, thyroid gland, parathyroid gland, salivary gland",
#               "skin, oral mucosa, adipose tissue",
#               "skeletal muscle, smooth muscle, heart muscle, cartilage, soft tissue",
#               "cerebellum, cerebral cortex, hippocampus, caudate, retina",
#               "lymph node, spleen, tonsil, bone marrow")
# )

# Creating dataframe with cellular locations
cell_structure_categories <- tribble(
  ~Category,               ~Term,
  "Cytoplasm and Associated Structures", "Cytosol",
  "Cytoplasm and Associated Structures", "Cytoplasmic bodies",
  "Cytoplasm and Associated Structures", "Lipid droplets",
  "Cytoplasm and Associated Structures", "Vesicles",
  "Cytoplasm and Associated Structures", "Aggresome",
  
  "Cytoskeleton", "Actin filaments",
  "Cytoskeleton", "Microtubules",
  "Cytoskeleton", "Intermediate filaments",
  "Cytoskeleton", "Mitotic spindle",
  "Cytoskeleton", "Mitotic chromosome",
  "Cytoskeleton", "Cytokinetic bridge",
  "Cytoskeleton", "Midbody",
  "Cytoskeleton", "Midbody ring",
  "Cytoskeleton", "Microtubule ends",
  "Cytoskeleton", "Rods & Rings",
  
  "Organelles", "Mitochondria",
  "Organelles", "Endoplasmic reticulum",
  "Organelles", "Golgi apparatus",
  "Organelles", "Peroxisomes",
  "Organelles", "Endosomes",
  "Organelles", "Lysosomes",
  "Organelles", "Centrosome",
  "Organelles", "Centriolar satellite",
  
  "Nucleus and Nuclear Structures", "Nucleoplasm",
  "Nucleus and Nuclear Structures", "Nuclear membrane",
  "Nucleus and Nuclear Structures", "Nucleoli",
  "Nucleus and Nuclear Structures", "Nucleoli fibrillar center",
  "Nucleus and Nuclear Structures", "Nucleoli rim",
  "Nucleus and Nuclear Structures", "Nuclear bodies",
  "Nucleus and Nuclear Structures", "Nuclear speckles",
  "Nucleus and Nuclear Structures", "Kinetochore",
  
  "Membrane-associated Structures", "Plasma membrane",
  "Membrane-associated Structures", "Cell Junctions",
  "Membrane-associated Structures", "Focal adhesion sites"
)



# Creating the dataframe using tribble
anatomical_categories <- tribble(
  ~Category,                    ~Term,
  # Digestive System
  "Digestive System",                "duodenum",
  "Digestive System",                "esophagus",
  "Digestive System",                "stomach",
  "Digestive System",                "small intestine",
  "Digestive System",                "colon",
  "Digestive System",                "rectum",
  "Digestive System",                "pancreas",
  "Digestive System",                "liver",
  "Digestive System",                "gallbladder",
  "Digestive System",                "appendix",
  "Digestive System",                "salivary gland",
  
  # Reproductive System - Female
  "Reproductive System - Female",    "breast",
  "Reproductive System - Female",    "cervix",
  "Reproductive System - Female",    "endometrium",
  "Reproductive System - Female",    "fallopian tube",
  "Reproductive System - Female",    "ovary",
  "Reproductive System - Female",    "placenta",
  "Reproductive System - Female",    "vagina",
  "Reproductive System - Female",    "lactating breast",
  
  # Reproductive System - Male
  "Reproductive System - Male",      "testis",
  "Reproductive System - Male",      "prostate",
  "Reproductive System - Male",      "epididymis",
  "Reproductive System - Male",      "seminal vesicle",
  
  # Respiratory System
  "Respiratory System",              "bronchus",
  "Respiratory System",              "lung",
  "Respiratory System",              "nasopharynx",
  
  # Nervous System - Central
  "Nervous System - Central",        "cerebral cortex",
  "Nervous System - Central",        "hippocampus",
  "Nervous System - Central",        "cerebellum",
  "Nervous System - Central",        "caudate",
  "Nervous System - Central",        "substantia nigra",
  "Nervous System - Central",        "dorsal raphe",
  "Nervous System - Central",        "hypothalamus",
  "Nervous System - Central",        "choroid plexus",
  
  # Nervous System - Peripheral
  "Nervous System - Peripheral",     "eye",
  "Nervous System - Peripheral",     "retina",
  
  # Endocrine System
  "Endocrine System",                "adrenal gland",
  "Endocrine System",                "thyroid gland",
  "Endocrine System",                "parathyroid gland",
  "Endocrine System",                "pituitary gland",
  
  # Immune and Lymphatic System
  "Immune and Lymphatic System",     "lymph node",
  "Immune and Lymphatic System",     "thymus",
  "Immune and Lymphatic System",     "spleen",
  "Immune and Lymphatic System",     "tonsil",
  "Immune and Lymphatic System",     "bone marrow",
  
  # Musculoskeletal System
  "Musculoskeletal System",          "skeletal muscle",
  "Musculoskeletal System",          "smooth muscle",
  "Musculoskeletal System",          "heart muscle",
  "Musculoskeletal System",          "cartilage",
  "Musculoskeletal System",          "soft tissue",
  "Musculoskeletal System",          "sole of foot",
  
  # Urinary System
  "Urinary System",                  "kidney",
  "Urinary System",                  "urinary bladder",
  
  # Integumentary System
  "Integumentary System",            "skin",
  "Integumentary System",            "hair",
  "Integumentary System",            "oral mucosa",
  "Integumentary System",            "adipose tissue",
)

