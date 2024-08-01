# Network_bio

This project analyzes the effect of ethoprophos on Daphnia magna gene expression levels collected at 3 different doses on 9 timepoints.

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
