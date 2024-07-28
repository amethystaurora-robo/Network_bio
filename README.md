# Network_bio

This project analyzes the effect of ethoprophos on Daphnia magna gene expression levels collected at 3 different doses on 9 timepoints.

Each code file contains a header which details its use.
The DynGENIE3 file requires a wrapper doc, available from https://github.com/vahuynh/dynGENIE3

DynGENIE3 is run in the following order:
Pre-processing_visualization.ipynb -> GRN_pre-processing.ipynb -> dyngenie_trial.R -> DynGENIE3_analysis.R -> process_results.R

Parameter tuning can be run on DynGENIE3 at any point after Pre-processing_visualization.ipynb. The parameter tuning file is parameter_tuning.R.

To annotate files and ensure robustness of DynGENIE3 analysis, Limma and GSEA are run. Limma is run using Galaxy software. GSEA is run using the online Gestalt tool. The files which pre-process the data for these tools are limma_pre-processing.ipynb and gestalt_pre-processing.ipynb.
