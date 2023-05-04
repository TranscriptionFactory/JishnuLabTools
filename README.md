# JishnuLabTools

## See the (properly formatted) HTML vignette
https://htmlpreview.github.io/?https://raw.githubusercontent.com/TranscriptionFactory/JishnuLabTools/main/running_ER_on_cluster.html?token=GHSAT0AAAAAACBC6C6R2H7JNAJYE5I3W5DMZCTXPOA


### Mostly tools for running Essential Regression (ER) and Significant Latent Factor Interaction Discovery and Exploration (SLIDE) 

Read the publications: 

1. **Essential Regression**
  - https://www.cell.com/patterns/pdfExtended/S2666-3899(22)00053-8
  - Github: See SLIDEpre repo below

2. **SLIDE**
  - (preprint) https://www.biorxiv.org/content/10.1101/2022.11.25.518001v1.full.pdf
  - Github: https://github.com/jishnu-lab/SLIDEpre



### Heuristics for running ER 
1. There are two important tuning variables: delta and lambda
2. Delta is more important and controls the number of latent factors
4. Lambda is less important and controls how latent factors are related

There are example yamls provided for regression and classification:
```
JishnuLabTools::regression_params
JishnuLabTools::classification_params
```
