### 🧬 Computational Modelling Pipeline for Omics Data  
This complete bioinformatics pipeline integrates traditional statistical analysis with cutting-edge graph machine learning to analyse transcriptomic and metabolomic samples of Daphnia magna. The organism in this pipeline has been treated with ethoprophos, an organophosphate pesticide, but the pipeline could be applied to any transcriptomic and metabolomic sample datasets. This allows for a network analysis of pathways affected under ethoprophos dosage.
<p align="left">
  <a href="https://public.tableau.com/app/profile/amethyst.eicher/vizzes" target="_blank">
    <img src="https://img.shields.io/badge/See%20Vizzes-766090?style=for-the-badge&logo=tableau&logoColor=white"/>
  </a>
</p>

### More Details Below
Each code file contains a header which details its use.
The DynGENIE3 file requires a wrapper doc, available from https://github.com/vahuynh/dynGENIE3

Pre-processed transcriptomic data is input into DynGENIE3, processed and annotated with Gestalt.
Pre-processing_visualization.ipynb -> GRN_pre-processing.ipynb -> dyngenie_trial.R -> DynGENIE3_analysis.R -> process_results.R -> gestalt_pre-processing.ipynb -> GESTALT using online software -> gsea_processing.ipynb -> network_processing.ipynb.

Pre-processed transcriptomic data and the csv output from network_processing.ipynb is input into:
WGCNA.r 
The csv output from this file can be used to plot the network in Cytoscape using module colors output from WGCNA. 

Parameter tuning can be run on DynGENIE3 at any point after Pre-processing_visualization.ipynb. The parameter tuning file is parameter_tuning.R.

DeSeq2 is run using raw transcriptomic data and files output from DeSeq2 in the following order:
rna_preproc.r -> rna_deseq.r -> gestalt_pre-processing.ipynb -> GESTALT using online software -> gsea_processing.ipynb 
The csv output is visualized using Tableau to see pathways significantly enriched over time and compare DEGs between DynGENIE3 and DeSeq2. 
