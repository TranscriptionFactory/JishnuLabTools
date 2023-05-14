# install.packages("WebGestaltR")
library(WebGestaltR)

# check for updates
# devtools::install_github("TranscriptionFactory/JishnuLabTools")
library(JishnuLabTools)
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

  # load the sig genes file
  sig_genes = JishnuLabTools:::load_output_from_runs(runs = er_input$out_path,
                                                     filename = "sig_genes")
}
