# Network_bio

This project analyzes the effect of ethoprophos on Daphnia magna gene expression levels collected at 3 different doses on 9 timepoints.

Each code file contains a header which details its use.
The DynGENIE3 file requires a wrapper doc, available from https://github.com/vahuynh/dynGENIE3

DynGENIE3 is run in the following order:
Pre-processing_visualization.ipynb -> GRN_pre-processing.ipynb -> dyngenie_trial.R -> DynGENIE3_analysis.R -> process_results.R

Parameter tuning can be run on DynGENIE3 at any point after Pre-processing_visualization.ipynb. The parameter tuning file is parameter_tuning.R.

DeSeq2 is run using raw transcriptomic data and files output from DeSeq2 and DynGENIE3 are processed for Gestalt.
rna_preproc.r -> rna_deseq.r -> gestalt_pre-processing.ipynb -> GESTALT using online software -> gsea_processing.ipynb
