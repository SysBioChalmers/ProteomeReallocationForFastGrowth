# ProteomeReallocationForFastGrowth
This repository contains a collection of supplementary data and scripts for analyzing and visualizing the data pertaining to the study:

_Proteome re-allocation from amino acid biosynthesis to ribosomes enables yeast to grow faster in rich media._

In this study, _Saccharomyces cerevisiae_ strain CEN.PK113-7D was cultivated in bioreactors with glucose as carbon source with or without supplementation of amino acids, both aerobically and anaerobically.

The main script **Organize_filter_and_quantify_proteomics.R** is used to normalize the TMT data into mass-percentages. The result of running this script is 6 data files named glu_proteomics_\*.csv that are used throughout the analysis. Input data for the script are **TMT_aerobic_cultures.csv**, **TMT_anaerobic_cultures.csv** and **ref2824_iBAQ_unique_peptides.csv**. For the iBAQ-file, fmol/ug amounts and Uniprot accessions were obtained from the proteomics analysis and matched to gene names and molecular weights obtaine from Uniprot. All original proteomics data can be obtained from the [PRIDE database](https://www.ebi.ac.uk/pride/) using accession numbers PXD018361 and PXD012803.

  ### Required software
  - R studio (version 1.0.136 or later)
