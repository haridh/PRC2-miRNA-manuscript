# PRC2-miRNA-manuscript
```
norm_counts_matrix_VP55_miRNA.R
```
Used for Fig. 2B. This script finds differentially expressed miRNAs using DESeq2 and makes a volcano plot. log2Pvalue vs log2foldchange.

```
norm_counts_matrix.R
norm_counts_matrix_VP55.R
```
These scripts identify differentially expressed mRNAs using DESeq2.
```
DESeq2_EZH2KO_WT_enrich_padj005_norm_counts_logfold_above_below_1_plus1_ratios.R
```
Heatmap for differentially expressed genes between EZH2-/- and WT. Fig. 3A
```
density_plot.R
```
Density plot fig 3C. Moving window average of ChIP-seq enrichment scores for EZH2 regulated genes.

```
EZH2_KO_WT_VP55_intersections.R
peak_side_analysis.R
```
Heatmap for miRNA repressed genes bound and regulated by EZH2. Fig. 3D and 4B.
```
scatter_ggplot_genes_analysis2.R
```
Scatter plot fig. 3F.
```
box_plot.R
```
Box plot for Fig.4C.

```
1000_shuffletest.R
miRNA_analysis_miR_Rshuffle.R
miRNA_analysis_target_Rshuffle.R
Interaction_map_final_gene_set_miR.R
```
Shuffle analysis. To determine background number of interactions for any 14 miRNA and 116 protein-coding genes. "Interaction_map_final_gene_set_miR.R" makes the plot in fig 4F.

```
scatter_ggplot_miRNA_analysis.R
```
Supplemental fig 1A. Scatter plot for miRNA enrichment in AGO2 iCLIP and input with miRNA labels.

```
Volcano_plot.R
```
Volcano plot for Supplemental fig. 3B.

```
Avg_plot_conf_intervals.R
```
ChIP-seq average plots for Supplemental fig. 6B

```
box_plot_ko_VP_in_set.R
```
Supplemental fig 8. 

```
box_plot_non_feed_for.R
```
Supplemental fig 9. 

```
box_plot_20_80_comp.R
Avg_plot_20_80_chip.R
```
Box plot and avg plot in Supplemental fig 11.
