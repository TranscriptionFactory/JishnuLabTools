# install.packages("WebGestaltR")
library(WebGestaltR)
library(tidyverse)
# check for updates
# devtools::install_github("TranscriptionFactory/JishnuLabTools")
library(JishnuLabTools)

# all organisms
organisms = listOrganism()

# picking mouse
organism = organisms[10]

# enrichment databases for gene names
gene_sets = (listGeneSet(organism = organism) %>%
  filter(idType == 'entrezgene') )$name

# all reference sets to use
reference_sets = listReferenceSet(organism = organism)

# genome protein coding
ref_set = reference_sets[25]

# all ID types (we're just using genename)
id_types = listIdType(organism = organism)

# using gene symbol
id = id_types[42]


run_params = list(
  organism = organism,
  gene_sets = gene_sets,
  ref_set = ref_set,
  id = id
)


get_enrichment_from_sig_genes = function(data_path, run_params,
                                    yaml_path = NULL) {

  # go to yaml out path and get the sig_genes file that has all the significant
  # latent factors
  if (!is.null(yaml_path) && file.exists(yaml_path)) {
    data_path = yaml::yaml.load_file(yaml_path)$out_path
  }

  # filename just needs to be a string specific to the output files you want
  sig_genes = JishnuLabTools:::load_output_from_runs(runs = data_path,
                                                     filename = "sig_genes")

  # get all the genes from all latent factors
  all_sig_genes = unname(unlist(lapply(unlist(sig_genes, recursive = F), function(x) x$gene)))

  # convert sig genes to title (first letter capitalized, rest lowercase)
  all_sig_genes = stringr::str_to_title(all_sig_genes)

  res <- WebGestaltR::WebGestaltR(enrichMethod = "ORA",
                                 organism = run_params$organism,
                                 interestGene = all_sig_genes,
                                 interestGeneType = run_params$id,
                                 enrichDatabase = run_params$gene_sets,
                                 referenceSet = run_params$ref_set,
                                 networkConstructionMethod = "Network_Retrieval_Prioritization",
                                 isOutput = TRUE,
                                 outputDirectory = data_path)
  # save RDS in the project folder
  data_path_dirs = list.dirs(data_path, full.names = T, recursive = F)

  # the projects are numbered in increasing order by default so the max # is the most
  # recent we just generated
  project_dir = max(data_path_dirs[stringr::str_which(data_path_dirs, pattern = "Project")])

  saveRDS(res, paste0(project_dir, "/results.RDS"))

  return(res)
}
