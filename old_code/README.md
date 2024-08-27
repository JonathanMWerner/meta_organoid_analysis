# meta_organoid_analysis
Repository for all the code used in the analysis of Werner and Gillis 2023

File descriptions

adult_brain_coexpression.Rmd - Used to generate the adult co-expression network and compare the preservation of co-expression of organoids between adult and fetal data. Contains graphs for Supplemental Figure 6 C - E

arlotta_organoid_atlas_analysis.Rmd - Normalization and QC for the organoid temporal atlas from the Arlotta lab

cell_type_specific_pres_coexp.Rmd - Computes the preserved co-expression of fine resolution cell-type markers in organoid data. Contains graphs for Figure 4 D and E

conservation_coexpression_fetal_organoids.Rmd - Computes the preservation of fetal co-expression in organoid co-expression networks. Contains graphs for Figure 2A, Figure 4 C and F, Supp. Fig. 6 E, and Supp. Fig. 7, 

dataset_download_exploration.Rmd - Collates metadata for the organoid datasets, performs normalization and QC for all organoid datasets and several annotated and all unannotated fetal datasets.

downsample_fetal_global_consCoexp.Rmd - Performs downsampling of fetal data to generate data for Supp. Fig. 6D

exploring_organoid_with_fetal_markers.Rmd - Collates Principle Components across the organoid datasets and explores the prescence of the fetal MetaMarkers. Contains graphs for Supp. Fig. 3C and D

fetal_coexpression.Rmd - Generates co-expression networks for the annotated and unannotated fetal data (excluding the first trimester Linnarsson dataset) and performs EGAD analysis for all fetal data.

fetal_conservation_coexpression_fetal.Rmd - Computes the preservation of fetal co-expression across fetal datasets. Contains graphs for Figure 1F, Figure 4B and G, Figure 6B, Supp. Fig. 5, Supp. Fig. 6B

fetal_data_with_annotations.Rmd - Collates the broad class level cell-type annotations for the annotated fetal datasets, the annotations used in generating MetaMarkers. Generates dataset specific markers for the annotated fetal datasets.

fetal_metaMarkers_cross_validation.Rmd - Performs the leave-one-out cross-validation of the fetal MetaMarkers. Contains graphs for Figure 2C-E, Supp. Figs. 1 - 2.

linnarsson_first_trimester_dataset.Rmd - Collates the broad class level cell-type annotations for the annotated first trimester fetal dataset from the Linnarsson lab. Performs normalization, QC, computes dataset specific cell-type markers, and computes the co-expression networks for this dataset.

organoid_coexpression.Rmd - Computes co-expression networks and performs EGAD analysis for all organoid datasets. Contains graphs for Figure 3B - E, and Supp. Fig. 4.

organoid_time_series_analysis.Rmd - Quantifies both predicting organoid cell-types and the preservation of co-expression of individual fetal dataset cell-type markers across organoid timepoints within the temporal organoid Arlotta atlas. Contains graphs for Figure 1B, Figure 5, and Supp. Fig. 3B  

process_fetal_areal_dataset.Rmd - Collates cell-type annotations and performs normalization and QC for the annotated fetal arealization datasets.

python_coexpression.py - python script for generating rank standardized co-expression networks

r_package_examples.Rmd - Compares the timing of computing co-expression networks and the preservation of co-expression across machines. Contains graphs for Figre 6A.
