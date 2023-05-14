# install.packages("WebGestaltR")
library(WebGestaltR)

# all organisms
organisms = listOrganism()

# picking mouse
organism = organisms[10]

gene_sets = listGeneSet(organism = organism)
reference_sets = listReferenceSet(organism = organism)
id_types = listIdType(organism = organism)




get_enrichment_from_yaml = function(yaml_path, organism = "hsapiens") {

  # go to yaml out path and get the sig_genes file that has all the significant
  # latent factors
  if (file.exists(yaml_path)) {
    er_input = yaml::yaml.load_file(yaml_path)
  }

}
